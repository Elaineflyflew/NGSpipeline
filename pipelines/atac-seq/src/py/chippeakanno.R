#!/usr/bin/env Rscript
################################################################################
# ⑫b ChIPpeakAnno —— 第二条注释路线
#
# 为什么要两条：
#   HOMER 用的是它自己整理的 hg38 注释（based on refGene），
#   ChIPpeakAnno 这里用的是 (1) 预制的 TSS.human.GRCh38 和 (2) UCSC knownGene。
#   两条路线给出的 promoter 占比接近，说明 peak 的位置没被某套注释带偏；
#   差得远就要回头看 ⑪ 的 TSS 图和 ⑦ 的片段分布。
#
# 用法：
#   Rscript chippeakanno.R --peaks x.narrowPeak --outdir 10_annot/chippeakanno --sample S1
################################################################################

suppressPackageStartupMessages({
  library(ChIPpeakAnno)
  library(GenomicFeatures)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1]
}

peakfile <- get_arg("--peaks")
outdir   <- get_arg("--outdir")
sample   <- get_arg("--sample", "sample")

if (is.null(peakfile) || is.null(outdir)) {
  stop("用法: Rscript chippeakanno.R --peaks <narrowPeak> --outdir <dir> --sample <name>")
}
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# narrowPeak 用 "narrowPeak" 解析，MACS2 的 *_peaks.xls 用 "MACS2"
fmt <- if (grepl("\\.narrowPeak$", peakfile)) "narrowPeak" else "MACS2"
cat(sprintf("[chippeakanno] 读 %s（format=%s）\n", peakfile, fmt))

peaks <- toGRanges(peakfile, format = fmt)
cat(sprintf("[chippeakanno] %d 个 peak\n", length(peaks)))

# 把注释结果写成「类别 + 个数 + 占比」三列，和 HOMER 的 annotation_distribution.tsv 同构
write_dist <- function(x, path) {
  tb <- table(x$insideFeature)
  df <- data.frame(feature = names(tb),
                   count = as.integer(tb),
                   pct = round(100 * as.integer(tb) / sum(tb), 2))
  write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(df)
}

# ---------------------------------------------------------------- 路线 1：预制 TSS
r1 <- tryCatch({
  data(TSS.human.GRCh38)
  anno <- annotatePeakInBatch(peaks, AnnotationData = TSS.human.GRCh38)
  anno <- addGeneIDs(annotatedPeak = anno, orgAnn = "org.Hs.eg.db",
                     IDs2Add = "symbol")
  write.table(as.data.frame(unname(anno)),
              file.path(outdir, paste0(sample, ".chippeakanno.tss.txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  pdf(file.path(outdir, paste0(sample, ".chippeakanno.tss.pie.pdf")), width = 7, height = 6)
  pie1(table(anno$insideFeature), main = "TSS.human.GRCh38")
  dev.off()
  write_dist(anno, file.path(outdir, paste0(sample, ".chippeakanno.tss.dist.tsv")))
}, error = function(e) {
  message("[chippeakanno] TSS.human.GRCh38 路线失败: ", conditionMessage(e))
  NULL
})

# ---------------------------------------------------------------- 路线 2：knownGene
r2 <- tryCatch({
  genes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
  anno <- annotatePeakInBatch(peaks, AnnotationData = genes)
  anno <- addGeneIDs(annotatedPeak = anno, orgAnn = "org.Hs.eg.db",
                     feature_id_type = "entrez_id", IDs2Add = "symbol")
  write.table(as.data.frame(unname(anno)),
              file.path(outdir, paste0(sample, ".chippeakanno.gene.txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  pdf(file.path(outdir, paste0(sample, ".chippeakanno.gene.pie.pdf")), width = 7, height = 6)
  pie1(table(anno$insideFeature), main = "UCSC hg38 knownGene")
  dev.off()
  write_dist(anno, file.path(outdir, paste0(sample, ".chippeakanno.gene.dist.tsv")))
}, error = function(e) {
  message("[chippeakanno] knownGene 路线失败: ", conditionMessage(e))
  NULL
})

# ---------------------------------------------------------------- 指标
metrics <- list()
if (!is.null(r1)) {
  metrics[["peaks_annotated"]] <- sum(r1$count)
  p <- r1$pct[r1$feature == "inside"]
  if (length(p) == 1) metrics[["inside_feature_pct"]] <- p
}
if (!is.null(r2)) {
  metrics[["peaks_annotated_gene"]] <- sum(r2$count)
}

con <- file(file.path(outdir, paste0(sample, ".chippeakanno_metrics.tsv")), "w")
writeLines("key\tvalue", con)
for (k in names(metrics)) writeLines(paste(k, metrics[[k]], sep = "\t"), con)
close(con)

cat(sprintf("[chippeakanno] 完成 → %s/\n", outdir))
