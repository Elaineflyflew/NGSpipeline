#!/usr/bin/env Rscript
################################################################################
# ⑭ ATACseqQC —— 用 Bioconductor 的实现把前面几步的关键结论再验一遍
#
# 为什么用**未偏移**的 srt.bam：
#   TSSE 和片段分布都基于"真实片段长度"。alignmentSieve --ATACshift 把每条 read
#   的坐标挪了 +4/-5，TLEN 字段却不会跟着改，两者会自相矛盾。
#   ATACseqQC 官方文档也要求输入未 shift 的 BAM。
#
# 这一步跑的东西（按顺序）：
#   1. bamQC            总 reads / 重复率 / 线粒体率 / 配对率 / MAPQ
#   2. estimateLibComplexity  文库复杂度（去重后还剩多少）
#   3. fragSizeDist     片段长度分布 —— 与 ⑦ 同口径
#   4. shiftGAlignmentsList   Tn5 偏移（ATACseqQC 自己算一遍）
#   5. PTscore / NFRscore     promoter/transcript 与 NFR 打分
#   6. TSSEscore        TSS 富集 —— 与 ⑪ 的 deepTools 结果对照
#   7. splitGAlignmentsByCut  NFR / 单 / 二 / 三核小体分档 + 累积占比
#   8. enrichedFragments + heatmap/distribution   TSS 两侧的核小体排布
#   9. factorFootprints  CTCF 足迹（可选，dbus 环境下 logo 可能画不出来）
#
# 设计原则：每一步都包在 tryCatch 里。⑭ 是交叉验证，任何一步失败都不该让
# 整个流程停下来——能出几个数就出几个数，最后在 TSV 里记下哪些步没跑成。
#
# 用法：
#   Rscript atacqc.R --bam srt.bam --outdir 12_atacqc --sample S1 \
#                    [--region ""] [--demo 1] [--footprint 1] [--cores 8]
################################################################################

suppressPackageStartupMessages({
  library(ATACseqQC)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(Rsamtools)
  library(GenomicRanges)
  library(GenomicAlignments)
  # estLibSize / reCenterPeaks / featureAlignedHeatmap / featureAlignedDistribution
  # 这几个看着像 ATACseqQC 的，其实都在 ChIPpeakAnno 里——ATACseqQC 的
  # vignette 是靠着同时 attach ChIPpeakAnno 才能跑通的，单独 library(ATACseqQC)
  # 会在第 8 步报 "could not find function"。
  library(ChIPpeakAnno)
})

# BSgenome 是第 7 步（splitGAlignmentsByCut 要 genome = Hsapiens）和第 9 步
# （factorFootprints）的硬需求。为什么单独软加载：服务器上常同时装着两套 R
# （系统 R 4.6 与 conda R 4.4），Bioconductor 包装在两套里的 libPaths 不同，
# 很容易出现"ATACseqQC 有、BSgenome 没有"。一个可选步骤不该把整步 ⑭ 带崩，
# 所以这里软探测，缺了就跳过 7/9 两步并记日志。
have_bsgenome <- suppressWarnings(suppressPackageStartupMessages(
  requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)))

# Rscript 把 stdout 重定向到文件时是块缓冲：不主动 flush 的话，日志要等这一大块
# 计算结束才整段冒出来，看上去像"卡死了"。每一步都是分钟级，必须逐行刷。
say <- function(...) { cat(...); flush(stdout()) }

if (!have_bsgenome) {
  say("[atacqc] 没装 BSgenome.Hsapiens.UCSC.hg38，第 9 步 CTCF 足迹会跳过\n")
} else {
  # 第 7 步要传 genome = Hsapiens，这个对象得先 attach 才在搜索路径上
  suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))
}

# ---------------------------------------------------------------- 参数
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1]
}

bamfile <- get_arg("--bam")
outdir  <- get_arg("--outdir")
sample  <- get_arg("--sample", "sample")
region  <- get_arg("--region", "")        # "" = 全基因组；"chr1" = 只读一条染色体
demo_region <- get_arg("--demo-region", "chr1")  # 核小体分档/heatmap 只在一条染色体上做
demo    <- get_arg("--demo", "1") == "1"  # 核小体分档 / heatmap / 足迹
footprint <- get_arg("--footprint", "1") == "1"

if (is.null(bamfile) || is.null(outdir)) {
  stop("用法: Rscript atacqc.R --bam <srt.bam> --outdir <dir> --sample <name>")
}
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
outPath <- file.path(outdir, "splited")
dir.create(outPath, showWarnings = FALSE, recursive = TRUE)

say(sprintf("[atacqc] 输入 %s\n", bamfile))
say(sprintf("[atacqc] 区域 %s\n", if (nzchar(region)) region else "全基因组"))

# ---------------------------------------------------------------- 指标收集
metrics <- list()
# 指标统一走这里。ATACseqQC 各版本返回的字段类型不一致（numeric / list 都有，
# 实测 1.30 的 qc$MAPQ 就是 list），所以先 unlist 再转 numeric，取第一个非 NA
# 值——一个怪字段不该把整步带崩。
add <- function(k, v) {
  if (is.null(v) || length(v) == 0) return(invisible(NULL))
  v <- suppressWarnings(as.numeric(unlist(v)))
  v <- v[!is.na(v)]
  if (length(v) == 0) return(invisible(NULL))
  metrics[[k]] <<- v[1]
  invisible(NULL)
}
# median 的安全版：先 unlist；空向量返回 NA（add() 会跳过）
med1 <- function(x) {
  x <- suppressWarnings(as.numeric(unlist(x)))
  if (length(x) == 0) return(NA_real_)
  stats::median(x, na.rm = TRUE)
}
# 每一步包一层：失败只记日志，不中断
step <- function(name, expr) {
  t0 <- Sys.time()
  d0 <- grDevices::dev.cur()          # 进来时当前设备，用来收尾
  res <- tryCatch(expr, error = function(e) {
    message(sprintf("[atacqc] %s 失败: %s", name, conditionMessage(e)))
    NULL
  })
  # 这一步里开的 pdf 设备在这里统一关掉。
  # 为什么不用 on.exit：expr 是 promise，在 step() 里被强制求值时它的
  # on.exit 落不到 step 的帧上，设备就一直开着——表现是 PDF 一直是 0 字节，
  # 要等 R 退出才落盘；万一某步中途报错还会串到别的图上。
  while (grDevices::dev.cur() > 1 && grDevices::dev.cur() != d0) {
    grDevices::dev.off()
  }
  dt <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  if (is.null(res)) {
    say(sprintf("[atacqc] %-22s FAIL (%ss)\n", name, dt))
  } else {
    say(sprintf("[atacqc] %-22s ok   (%ss)\n", name, dt))
  }
  invisible(res)
}

pdf_open  <- function(f, w = 8, h = 5) pdf(f, width = w, height = h)
pdf_close <- function() if (dev.cur() > 1) dev.off()

# ---------------------------------------------------------------- 1. bamQC
# 25 M reads 的 BAM 上 bamQC 要 8 分钟左右（单线程读一遍算重复率/线粒体率/
# MAPQ 分布），所以别指望它快。
# 注意：这里读的是 ④ 过滤后的 srt.bam，chrM 在 ③ 的比对过滤里就被剔掉了，
# 所以 mitochondria_rate 会算成 0 —— 真实的 chrM 比例看 ③ 的指标，
# ⑭ 这个 0 只是"过滤确实生效了"的旁证。
qc <- step("bamQC", bamQC(bamfile))
if (!is.null(qc)) {
  # 取指标本身也包一层：字段类型在不同版本里会变，出错只丢指标不丢流程
  step("bamQC-metrics", {
    add("total_qnames",       med1(qc$totalQNAMEs))
    add("duplicate_rate",     round(med1(qc$duplicateRate), 4))
    add("mitochondria_rate",  round(med1(qc$mitochondriaRate), 4))
    add("proper_pair_rate",   round(med1(qc$properPairRate), 4))
    add("unmapped_rate",      round(med1(qc$unmappedRate), 4))
    add("non_redundant_frac", round(med1(qc$nonRedundantFraction), 4))
    # qc$MAPQ 不是向量，是 data.frame(Var1 = MAPQ 值, Freq = 条数)。
    # 直接 unlist 会先取到 Freq 列，算出来是"某个 MAPQ 分箱的条数"（实测 1956），
    # 必须按条数加权求中位数。
    mq <- qc$MAPQ
    if (is.data.frame(mq) && all(c("Var1", "Freq") %in% names(mq))) {
      v <- suppressWarnings(as.numeric(as.character(mq$Var1)))
      w <- suppressWarnings(as.numeric(mq$Freq))
      ok <- !is.na(v) & !is.na(w) & w > 0
      if (any(ok)) {
        o <- order(v[ok])
        v <- v[ok][o]
        w <- w[ok][o]
        add("median_mapq", v[which(cumsum(w) >= sum(w) / 2)[1]])
      }
    } else {
      add("median_mapq", med1(mq))
    }
    # idxstats 里有每条染色体的比对 reads 数，单独存一份方便查
    if (!is.null(qc$idxstats)) {
      try(write.table(qc$idxstats,
                      file.path(outdir, paste0(sample, ".idxstats.tsv")),
                      sep = "\t", quote = FALSE, row.names = FALSE,
                      col.names = FALSE), silent = TRUE)
    }
    TRUE
  })
}

# ---------------------------------------------------------------- 2. 文库复杂度
step("libComplexity", {
  pdf_open(file.path(outdir, paste0(sample, ".libComplexity.pdf")))
  estimateLibComplexity(readsDupFreq(bamfile))
  TRUE
})

# ---------------------------------------------------------------- 3. 片段分布
step("fragSizeDist", {
  pdf_open(file.path(outdir, paste0(sample, ".fragSizeDist.pdf")))
  fragSizeDist(bamfile, gsub(".bam$", "", basename(bamfile)))
  TRUE
})

# ---------------------------------------------------------------- 4. 读入 + Tn5 偏移
possibleTag <- list(
  "integer" = c("AM", "AS", "CM", "CP", "FI", "H0", "H1", "H2", "HI", "IH",
                "MQ", "NH", "NM", "OP", "PQ", "SM", "TC", "UQ"),
  "character" = c("BC", "BQ", "BZ", "CB", "CC", "CO", "CQ", "CR", "CS", "CT",
                  "CY", "E2", "FS", "LB", "MC", "MD", "MI", "OA", "OC", "OQ",
                  "OX", "PG", "PT", "PU", "Q2", "QT", "QX", "R2", "RG", "RX",
                  "SA", "TS", "U2")
)
tags <- step("scanTags", {
  t0 <- scanBam(BamFile(bamfile, yieldSize = 100),
                param = ScanBamParam(tag = unlist(possibleTag)))[[1]]$tag
  names(t0)[lengths(t0) > 0]
})
if (is.null(tags)) tags <- character(0)

seqinfo_txdb <- seqinfo(TxDb.Hsapiens.UCSC.hg38.knownGene)
if (nzchar(region)) {
  if (!region %in% seqlevels(seqinfo_txdb)) {
    stop(sprintf("region %s 不在 TxDb 的染色体列表里", region))
  }
  which_gr <- as(seqinfo_txdb[region], "GRanges")
  gal <- step("readBamFile", readBamFile(bamfile, tag = tags, which = which_gr,
                                         asMates = TRUE, bigFile = FALSE))
} else {
  # 不给 which 就是全基因组。Rscript 里 interactive() 是 FALSE，
  # 所以 readBamFile 不会弹 "bigFile=TRUE?" 那个问题。
  gal <- step("readBamFile", readBamFile(bamfile, tag = tags,
                                         asMates = TRUE, bigFile = FALSE))
}

gal1 <- NULL
# 两个"版本坑"都在这一个调用上：
#  1) outbam 从 1.38 起是**必填**参数（没有默认值），不传直接报错；
#  2) 返回值从 GAlignmentsList 变成了**扁平的 GAlignments**（readBamFile 出来的
#     还是 list，shift 完就 unlist 了）。后面所有按 seqnames 取子集、按片段长度
#     分档的代码都必须按 GAlignments 写，不能再 endoapply。
# 产物我们不保留：足迹走 ④ 的 shift.bam（见第 9 步），这里写完就删，
# 25 M reads 的中间文件能到 5 GB。
shiftedBam <- file.path(outPath, "shifted.bam")
if (!is.null(gal)) {
  gal1 <- step("shiftGAlignments", shiftGAlignmentsList(gal, outbam = shiftedBam))
  for (f in c(shiftedBam, sub("[.]bam$", ".sam", shiftedBam),
              paste0(shiftedBam, ".bai"))) {
    if (file.exists(f)) unlink(f)
  }
}
if (!is.null(gal1)) {
  say(sprintf("[atacqc] shift 后对象类型 %s，%d 条 alignment\n",
              paste(class(gal1), collapse = "/"), length(gal1)))
}

txs <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)

# ---------------------------------------------------------------- 5. PT / NFR 打分
if (!is.null(gal1)) {
  step("PTscore", {
    pt <- PTscore(gal1, txs)
    add("pt_score_median", round(med1(pt$PT_score), 3))
    pdf_open(file.path(outdir, paste0(sample, ".PTscore.pdf")))
    plot(pt$log2meanCoverage, pt$PT_score,
         xlab = "log2 mean coverage", ylab = "Promoter vs Transcript")
    TRUE
  })

  step("NFRscore", {
    nfr <- NFRscore(gal1, txs)
    add("nfr_score_median", round(med1(nfr$NFR_score), 3))
    pdf_open(file.path(outdir, paste0(sample, ".NFRscore.pdf")))
    plot(nfr$log2meanCoverage, nfr$NFR_score,
         xlab = "log2 mean coverage", ylab = "Nucleosome Free Regions score",
         main = "NFRscore for 200bp flanking TSSs",
         xlim = c(-10, 0), ylim = c(-5, 5))
    TRUE
  })

  # ------------------------------------------------------------ 6. TSSE
  step("TSSEscore", {
    tsse <- TSSEscore(gal1, txs)
    add("tss_enrichment_score_r", round(med1(tsse$TSSEscore), 3))
    pdf_open(file.path(outdir, paste0(sample, ".TSSEscore.pdf")))
    plot(100 * (-9:10 - .5), tsse$values, type = "b",
         xlab = "distance to TSS", ylab = "aggregate TSS score")
    TRUE
  })
}

# ---------------------------------------------------------------- 7. 核小体分档
# 片段分档用的是 ATACseqQC 自己的阈值（NFR<100, mono 180-247, di 315-473, ...），
# 和 ⑦ 的 147/294 不完全一样——这正是"换个实现再算一遍"的意义：
# 两套阈值下 NFR 都是最大的那一档，才说明片段分布真的正常。
objs <- NULL
# 这一步要 genome（BSgenome）才能算 GC/保守性，没装就整段跳过——
# 分档本身能跑，但 ATACseqQC 的接口把 genome 设成必需参数，没法只给一半。
if (demo && have_bsgenome && !is.null(gal1)) {
  # 分档 + conservation 是逐 fragment 的操作，全基因组要跑很久。
  # 默认只取一条染色体做这套"演示"，TSSE 仍然用全基因组（那才是要交叉验证的数）。
  # 注意 gal1 是扁平的 GAlignments（见第 4 步），直接按 seqnames 取子集即可，
  # 不能用 endoapply —— 那是对 GAlignmentsList 的写法，在这里会报
  # "GAlignments objects don't support [[, $, as.list(), lapply()"。
  gal_demo <- gal1
  if (nzchar(demo_region)) {
    gal_demo <- gal1[as.character(seqnames(gal1)) %in% demo_region]
  }
  say(sprintf("[atacqc] 分档区域 %s，%d 条 alignment\n",
              if (nzchar(demo_region)) demo_region else "全基因组", length(gal_demo)))

  # conservation 是可选的（要 phastCons100way.UCSC.hg38 这个数据包，
  # 体积不小，很多集群没装）。没有就传 NULL，分档照样出。
  have_phast <- suppressWarnings(suppressPackageStartupMessages(
    requireNamespace("phastCons100way.UCSC.hg38", quietly = TRUE)))
  cons <- if (have_phast) {
    phastCons100way.UCSC.hg38::phastCons100way.UCSC.hg38
  } else {
    say("[atacqc] 没装 phastCons100way.UCSC.hg38，分档不做保守性注释\n")
    NULL
  }

  objs <- step("splitByCut", {
    splitGAlignmentsByCut(gal_demo, txs = txs, genome = Hsapiens, outPath = outPath,
                          conservation = cons)
  })

  bamfiles <- file.path(outPath,
                        c("NucleosomeFree.bam", "mononucleosome.bam",
                          "dinucleosome.bam", "trinucleosome.bam"))
  bamfiles <- bamfiles[file.exists(bamfiles)]

  # 分档占比（与 ⑦ 的 147/294 不是同一套阈值，正好用来交叉验证）
  if (!is.null(objs)) {
    cnt <- function(k) if (is.null(objs[[k]])) 0 else length(objs[[k]])
    tot <- sum(vapply(c("NucleosomeFree", "mononucleosome", "dinucleosome",
                        "trinucleosome"), cnt, numeric(1)))
    if (tot > 0) {
      add("nfr_fraction_r",  round(cnt("NucleosomeFree") / tot, 4))
      add("mono_fraction_r", round(cnt("mononucleosome") / tot, 4))
      add("poly_fraction_r", round((cnt("dinucleosome") + cnt("trinucleosome")) / tot, 4))
      add("split_fragments", tot)
    }
  }

  # ATACseqQC 1.38 把 cumulativePercentage() 删了（老版本有），所以这里自己画：
  # 横轴片段长度，纵轴"≤ 该长度的片段占本档的比例"。四条曲线的拐点位置就是
  # ⑦ 那张片段分布图的累积版——NFR 曲线最先冲到 1，mono 其次。
  step("cumulativePct", {
    bf <- bamfiles[grepl("NucleosomeFree|mononucleosome|dinucleosome|trinucleosome",
                         basename(bamfiles))]
    if (length(bf) == 0) stop("没有分档 BAM，跳过累积曲线")
    tlen <- lapply(bf, function(f) {
      abs(unlist(scanBam(f, param = ScanBamParam(what = "isize"))[[1]]$isize))
    })
    names(tlen) <- sub("[.]bam$", "", basename(bf))
    xs <- seq(0, 800, by = 5)
    cols <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")
    pdf_open(file.path(outdir, paste0(sample, ".cumulativePercentage.pdf")))
    plot(NA, xlim = c(0, 800), ylim = c(0, 1),
         xlab = "Fragment length (bp)", ylab = "Cumulative fraction")
    for (i in seq_along(tlen)) {
      v <- tlen[[i]]
      if (length(v) == 0) next
      lines(xs, vapply(xs, function(x) mean(v <= x), numeric(1)),
            col = cols[(i - 1) %% length(cols) + 1], lwd = 2)
    }
    legend("bottomright", legend = names(tlen),
           col = cols[seq_along(tlen)], lwd = 2, bty = "n", cex = 0.8)
    TRUE
  })

  # ------------------------------------------------------------ 8. TSS 两侧排布
  step("enrichedFragments", {
    TSS <- unique(promoters(txs, upstream = 0, downstream = 1))
    libSize <- estLibSize(bamfiles)
    NTILE <- 101; dws <- 1010; ups <- 1010
    keep <- intersect(c("NucleosomeFree", "mononucleosome", "dinucleosome",
                        "trinucleosome"), names(objs))
    sigs <- enrichedFragments(gal = objs[keep], TSS = TSS, librarySize = libSize,
                              seqlev = if (nzchar(demo_region)) demo_region else seqlevels(seqinfo_txdb),
                              TSS.filter = 0.5, n.tile = NTILE,
                              upstream = ups, downstream = dws)
    sigs.log2 <- lapply(sigs, function(.ele) log2(.ele + 1))

    pdf_open(file.path(outdir, paste0(sample, ".featureAlignedHeatmap.pdf")))
    featureAlignedHeatmap(sigs.log2, reCenterPeaks(TSS, width = ups + dws),
                          zeroAt = .5, n.tile = NTILE)
    dev.off()

    pdf_open(file.path(outdir, paste0(sample, ".featureAlignedDistribution.pdf")))
    out <- featureAlignedDistribution(sigs, reCenterPeaks(TSS, width = ups + dws),
                                      zeroAt = .5, n.tile = NTILE, type = "l",
                                      ylab = "Averaged coverage")
    # 归一到 0~1 再画，四条曲线才在同一个量纲上
    range01 <- function(x) (x - min(x)) / (max(x) - min(x))
    out <- apply(out, 2, range01)
    matplot(out, type = "l", xaxt = "n",
            xlab = "Position (bp)", ylab = "Fraction of signal")
    axis(1, at = seq(0, 100, by = 10) + 1,
         labels = c("-1K", seq(-800, 800, by = 200), "1K"), las = 2)
    abline(v = seq(0, 100, by = 10) + 1, lty = 2, col = "gray")
    dev.off()
    TRUE
  })
}

# ---------------------------------------------------------------- 9. CTCF 足迹
# 为什么不复用内存里的 gal1：ATACseqQC ≥1.34 的 factorFootprints 只收磁盘上的
# BAM（要能建索引），不收 GAlignmentsList。所以这里直接用 ④ 里
# alignmentSieve --ATACshift 已经产出的 shift.bam，不再另写一份。
# 另外 ≥1.34 把 bindingSites 改成了必填参数（老版本会自己扫全基因组），
# 于是这里从 ⑨ 的 narrowPeak 里取区间当候选位点——本来也更有意义：
# 只在自己 call 出来的开放区里找 CTCF 足迹。
footprint_bam <- get_arg("--footprint-bam", "")
peaks_file    <- get_arg("--peaks", "")

if (footprint) {
  miss <- c(
    if (!have_bsgenome) "BSgenome.Hsapiens.UCSC.hg38",
    if (!nzchar(footprint_bam) || !file.exists(footprint_bam)) "带索引的 shift.bam",
    if (!nzchar(peaks_file) || !file.exists(peaks_file)) "⑨ 的 narrowPeak"
  )
  if (length(miss) > 0) {
    say(sprintf("[atacqc] 第 9 步 CTCF 足迹跳过，缺: %s\n", paste(miss, collapse = ", ")))
  } else {
    step("factorFootprints", {
      suppressPackageStartupMessages(library(MotifDb))
      if (!file.exists(paste0(footprint_bam, ".bai"))) {
        say("[atacqc] 给 shift.bam 建索引…\n")
        Rsamtools::indexBam(footprint_bam)
      }
      # narrowPeak 无表头，前三列就是区间
      pk <- utils::read.delim(peaks_file, header = FALSE, comment.char = "#")
      sites <- GRanges(seqnames = pk[[1]],
                       ranges = IRanges(start = pk[[2]] + 1L, end = pk[[3]]))
      if (nzchar(demo_region)) {
        sites <- sites[as.character(seqnames(sites)) %in% demo_region]
      }
      # factorFootprints 要求 bindingSites 带 seqlengths（它拿这个算期望位点数），
      # 而从 narrowPeak 现搭的 GRanges 没有 —— 从 BSgenome 补上，否则报
      # "不是所有的 all(!is.na(seqlengths(bindingSites))) 都是 TRUE"。
      sites <- GenomeInfoDb::keepSeqlevels(
        sites, intersect(seqlevels(sites), seqlevels(Hsapiens)),
        pruning.mode = "coarse")
      seqinfo(sites) <- seqinfo(Hsapiens)[seqlevels(sites)]
      say(sprintf("[atacqc] 足迹候选位点 %d 个（%s）\n", length(sites),
                  if (nzchar(demo_region)) demo_region else "全基因组"))
      CTCF <- as.list(MotifDb::query(MotifDb, c("CTCF")))
      # 无图形环境下 dbus 会报错（logo 画不出来），但曲线本身是好的
      Sys.setenv(NO_AT_BRIDGE = "1", DBUS_SESSION_BUS_ADDRESS = "",
                 QT_QPA_PLATFORM = "offscreen")
      pdf_open(file.path(outdir, paste0(sample, ".factorFootprints.pdf")))
      factorFootprints(footprint_bam, pfm = CTCF[[1]], genome = Hsapiens,
                       min.score = "90%", bindingSites = sites,
                       # 默认 seqlev 是 chr1..chr22+X+Y，而 BAM 头里未必有 chrY
                       # （这个样本 chrY 没 reads，header 里可能就没有），
                       # 那样会报 "There is NA values in seqlengths ... chrY"。
                       # 所以显式限定成候选位点所在的染色体。
                       seqlev = if (nzchar(demo_region)) demo_region
                                else as.character(unique(seqnames(sites))),
                       upstream = 100, downstream = 100)
      TRUE
    })
  }
}

# ---------------------------------------------------------------- 写指标
out_tsv <- file.path(outdir, paste0(sample, ".atacqc_metrics.tsv"))
con <- file(out_tsv, "w")
writeLines("key\tvalue", con)
for (k in names(metrics)) writeLines(paste(k, metrics[[k]], sep = "\t"), con)
close(con)

say(sprintf("\n[atacqc] %d 个指标 → %s\n", length(metrics), out_tsv))
say(sprintf("[atacqc] 图 → %s/*.pdf\n", outdir))
