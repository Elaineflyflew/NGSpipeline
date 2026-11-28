#!/usr/bin/env bash
# =============================================================================
#  src/10_cli.sh —— 命令行接口
#
#  只做两件事：打印帮助（usage / list_stages）、把 argv 解析成全局变量
#  （parse_cli）。新增一个参数 = usage 里加一行 + parse_cli 的 case 加一个分支。
# =============================================================================

# 全局状态（由 parse_cli 填充）
R1=""; R2=""; SAMPLE=""
STAGES=""
FORCE=0; DRY_RUN=0; CHECK_ONLY=0; KEEP_TMP=0; LIST_STAGES=0
CONFIG_FILE=""

usage() {
cat <<EOF
${PROG} ${VERSION} —— ATAC-seq 双端 FASTQ 分析流水线

用法
  ${PROG} -1 R1.fastq.gz -2 R2.fastq.gz [选项]

输入
  -1, --r1 FILE          双端 R1（.fastq/.fq，可 .gz）
  -2, --r2 FILE          双端 R2
  -s, --sample NAME      样本名（默认从文件名推导）
  -o, --outdir DIR       输出目录（默认: ${OUTDIR}）
  -t, --threads N        线程数（默认: ${THREADS}）

参考数据（默认都指向项目内 reference/，deploy 之后不用改）
  -g, --genome FILE      基因组 fasta（默认: reference/genome/ 下自动探测）
      --bt2 PREFIX       bowtie2 索引前缀（默认: reference/genome/bt2/ 下自动探测）
      --gtf FILE         基因注释 GTF（默认: reference/annotation/ 下第一个）
      --blacklist FILE   ENCODE 黑名单 BED（默认: reference/blacklist/ 下第一个）
      --homer-genome STR HOMER 基因组名（默认 ${HOMER_GENOME}）
      --homer-motifs FILE HOMER known.motifs（默认 ${HOMER_MOTIFS}）

阈值
      --mapq N           比对质量下限（默认 ${MIN_MAPQ}）
      --nfr-max N        NFR 上界 bp（默认 ${TLEN_NFR_MAX}）
      --mono-max N       单核小体上界 bp（默认 ${TLEN_MONO_MAX}）
      --qvalue X         MACS2 q 值（默认 ${MACS2_QVALUE}）
      --topn N           motif 用的 top peak 数（默认 ${MOTIF_TOPN}）

流程控制
      --stages LIST      qc,trim,align,bam,frag,peaks,frip,coverage,annot,motif,atacqc,igv,report
      --list-stages      打印所有阶段及产物后退出
      --force            忽略完成标记，强制重跑
      --check            只检查工具和参考数据，不跑数据
      --dry-run          只打印命令，不执行
      --keep-tmp         保留 results/.tmp（排序列临时文件），默认跑完清掉
  -C, --config FILE      额外加载一个配置文件（覆盖 config/atac.conf）
  -h, --help             显示本帮助
      --version          显示版本

IGV
      --igv-locus STR    截图区域（默认 ${IGV_LOCUS}，即 MYC）
      --igv-width N      截图宽（默认 ${IGV_WIDTH}）
      --igv-height N     截图高（默认 ${IGV_HEIGHT}）

环境变量
  上面每个参数都有同名大写变量（GENOME / BT2 / THREADS / OUTDIR ...），
  可以写进 ~/.bashrc 或 env 文件。

示例
  # 一条命令跑完 15 步
  ${PROG} -1 raw/sub20M.SRR13849313_1.fastq.gz -2 raw/sub20M.SRR13849313_2.fastq.gz -t 20

  # 上服务器前先看会跑什么
  ${PROG} -1 ... -2 ... --dry-run

  # 只重跑 peak 之后的部分
  ${PROG} -1 ... -2 ... --stages peaks,frip,coverage,annot,motif,atacqc,igv,report --force
EOF
}

list_stages() {
# 三列布局：内容放最后一列，不靠空格对齐（中文宽度算 2，heredoc 里对不准）
cat <<'EOF'
阶段     步骤   内容 → 产物
qc       ①      FastQC 原始质检 → 01_fastqc/*.html
trim     ②      fastp 接头修剪 + 清洗后复检 → 02_fastp/*.clean_*.fq.gz
align    ③      bowtie2 --local --very-fast-local → 03_map/*.sam
bam      ④⑤⑥   MAPQ 过滤 + 去 chrM + 排序 → Tn5 偏移 → 重排 → 04_bam/*.srt.bam, *.srt2.bam
frag     ⑦      按 |TLEN| 分 NFR / 单核小体 / 多核小体 → 05_fragbin/*.bam + fragsize_dist.png
peaks    ⑧⑨    MACS2 ×3 路线 + 去 ENCODE 黑名单 → 06_peaks/*.filtered.narrowPeak
frip     ⑩      FRiP 双算法交叉验证 → 07_frip/<sample>/*.frip.tsv
coverage ⑪      bamCoverage CPM + TSS 富集 → 08_bigwig/*.bw, 09_tss/*
annot    ⑫      HOMER annotatePeaks 基因组分布 → 10_annot/*
motif    ⑬      HOMER findMotifsGenome -mknown → 11_motif/*
atacqc   ⑭      ATACseqQC 片段分布 / TSSE → 12_atacqc/*
igv      ⑮      MYC 位点截图（IGV 或 pyGenomeTracks）→ 13_igv/*.png
report   ★      指标汇总 + PASS/WARN 判读 + HTML 报告 → 99_report/*
EOF
}

# 把 --opt=val 规整成 --opt val，后面只写一套 case（bash 3.2 兼容）
parse_cli() {
  if [[ $# -gt 0 ]]; then
    local -a _new=()
    local _a
    for _a in "$@"; do
      if [[ "$_a" == --*=* ]]; then _new+=("${_a%%=*}" "${_a#*=}"); else _new+=("$_a"); fi
    done
    set -- "${_new[@]}"
  fi

  while [[ $# -gt 0 ]]; do
    case "$1" in
      -1|--r1)          R1="$2"; shift 2 ;;
      -2|--r2)          R2="$2"; shift 2 ;;
      -s|--sample)      SAMPLE="$2"; shift 2 ;;
      -o|--outdir)      OUTDIR="$2"; shift 2 ;;
      -t|--threads)     THREADS="$2"; shift 2 ;;
      -g|--genome)      GENOME="$2"; shift 2 ;;
      --bt2)            BT2="$2"; shift 2 ;;
      --gtf)            GTF="$2"; shift 2 ;;
      --blacklist)      BLACKLIST="$2"; shift 2 ;;
      --homer-genome)   HOMER_GENOME="$2"; shift 2 ;;
      --homer-motifs)   HOMER_MOTIFS="$2"; shift 2 ;;
      --mapq)           MIN_MAPQ="$2"; shift 2 ;;
      --nfr-max)        TLEN_NFR_MAX="$2"; shift 2 ;;
      --mono-max)       TLEN_MONO_MAX="$2"; shift 2 ;;
      --qvalue)         MACS2_QVALUE="$2"; shift 2 ;;
      --topn)           MOTIF_TOPN="$2"; shift 2 ;;
      --stages)         STAGES="$2"; shift 2 ;;
      --igv-locus)      IGV_LOCUS="$2"; shift 2 ;;
      --igv-width)      IGV_WIDTH="$2"; shift 2 ;;
      --igv-height)     IGV_HEIGHT="$2"; shift 2 ;;
      -C|--config)      CONFIG_FILE="$2"; shift 2 ;;
      --force)          FORCE=1; shift ;;
      --check)          CHECK_ONLY=1; shift ;;
      --dry-run)        DRY_RUN=1; shift ;;
      --keep-tmp)       KEEP_TMP=1; shift ;;
      --list-stages)    LIST_STAGES=1; shift ;;
      -h|--help)        usage; exit 0 ;;
      --version)        printf '%s %s\n' "$PROG" "$VERSION"; exit 0 ;;
      *) printf '未知参数: %s\n\n' "$1" >&2; usage >&2; exit 2 ;;
    esac
  done
}
