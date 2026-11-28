#!/usr/bin/env bash
# =============================================================================
#  src/10_cli.sh —— 命令行接口
#
#  只做两件事：打印帮助（usage）、把 argv 解析成全局变量（parse_cli）。
#  新增一个参数 = usage 里加一行 + parse_cli 的 case 里加一个分支。
# =============================================================================

R1=""; R2=""; BAM_IN=""; SAMPLE_OPT=""
STAGES=""
FORCE=0; DRY_RUN=0; CHECK_ONLY=0; KEEP_TMP=0; LIST_STAGES=0
CONFIG_FILE=""

usage() {
cat <<EOF
${PROG} ${VERSION} —— ATAC-seq 开放染色质分析流水线

用法
  ${PROG} -1 R1.fastq.gz -2 R2.fastq.gz [选项]     # 从原始 FASTQ 开始
  ${PROG} --bam sample.srt2.bam        [选项]     # 复用已偏移的 BAM，直接做 peak

输入
  -1, --r1 FILE          双端 R1（.fastq/.fq，可 .gz）
  -2, --r2 FILE          双端 R2
      --bam FILE         已偏移的 BAM（srt2），跳过 qc/trim/align/bam/shift

核心参数
  -g, --genome NAME      物种配置（默认 ${GENOME}，读 config/genomes/NAME.conf）
  -o, --outdir DIR       输出目录（默认 ${OUTDIR}）
  -s, --sample NAME      样本名（默认从文件名推导）
  -t, --threads N        线程数（默认 ${THREADS}）

参考数据（默认从 config/genomes/${GENOME}.conf 探测，一般不用改）
      --bt2-index PREFIX   bowtie2 索引前缀
      --fasta FILE         参考基因组 fasta
      --blacklist FILE     ENCODE blacklist bed
      --gtf FILE           基因注释 gtf
      --motif-db FILE      JASPAR meme 库
      --homer-motifs FILE  HOMER known.motifs

阈值
      --min-mapq N         比对质量下限（默认 ${MIN_MAPQ}）
      --keep-mt            保留线粒体 reads（默认丢弃）
      --nfr-max N          NFR 上限 bp（默认 ${NFR_MAX}）
      --mono-max N         单核小体上限 bp（默认 ${MONO_MAX}）
      --macs2-qvalue X     MACS2 q 值（默认 ${MACS2_QVALUE}）
      --toppeaks N         motif 分析取前 N 个 peak（默认 ${TOPPEAKS_N}）

流程控制
      --stages LIST        qc,trim,align,bam,shift,fraction,peak,metrics,bigwig,annotate,atacqc,report
      --list-stages        打印所有阶段及说明后退出
      --force              忽略完成标记，强制重跑
      --check              只检查工具和参考数据，不跑数据
      --dry-run            只打印命令，不执行
      --keep-tmp           保留中间文件
      --no-motif           跳过 ⑬ motif 富集
      --no-atacqc          跳过 ⑭ 片段分布 / TSSE
      --meme               额外跑 MEME motif（很慢，默认关）
      --tobias             额外跑 TOBIAS footprint
  -C, --config FILE        额外加载一个配置文件
  -h, --help               显示本帮助
      --version            显示版本

IGV
      --igv-region STR     截图区域（默认 ${IGV_REGION}，即 MYC）

示例
  # 全流程
  ${PROG} -1 raw/Case1_R1.fastq.gz -2 raw/Case1_R2.fastq.gz -t 20

  # 换物种
  ${PROG} -1 ... -2 ... -g mm10

  # 只重跑 peak calling
  ${PROG} --bam results/02_bamfiles/Case1.srt2.bam --stages peak,report --force
EOF
}

list_stages() {
cat <<EOF
可用阶段（按执行顺序）：

  qc        ① FastQC 原始数据质控
  trim      ② fastp 接头修剪 + 质量过滤
  align     ③ bowtie2 比对 → SAM
  bam       ④ samtools view -q / 去 chrM → sort → index
  shift     ⑤⑥ alignmentSieve --ATACshift + 重排建索引
  fraction  ⑦ 按 TLEN 分档：NFR / 单核小体 / 多核小体
  peak      ⑧⑨ MACS2 callpeak ×3 + bedtools 去 blacklist
  metrics   ⑩ FRiP 双算法交叉验证 + TSS 富集分数
  bigwig    ⑪ bamCoverage → bigWig，TSS 富集曲线
  annotate  ⑫⑬ HOMER peak 注释 + motif 富集
  atacqc    ⑭ 片段长度分布 + TSSE（内置实现）
  report    ⑮ 汇总指标 + 自动结论 + IGV 批处理
EOF
}

parse_cli() {
  while (( $# )); do
    case "$1" in
      -1|--r1)          R1="$2"; shift 2 ;;
      -2|--r2)          R2="$2"; shift 2 ;;
      --bam)            BAM_IN="$2"; shift 2 ;;
      -g|--genome)      load_genome "$2"; shift 2 ;;
      -o|--outdir)      OUTDIR="$2"; shift 2 ;;
      -s|--sample)      SAMPLE_OPT="$2"; shift 2 ;;
      -t|--threads)     THREADS="$2"; shift 2 ;;
      --bt2-index)      BT2_INDEX="$2"; shift 2 ;;
      --fasta)          FASTA="$2"; shift 2 ;;
      --blacklist)      BLACKLIST="$2"; shift 2 ;;
      --gtf)            GTF="$2"; shift 2 ;;
      --motif-db)       MOTIF_DB="$2"; shift 2 ;;
      --homer-motifs)   HOMER_MOTIFS="$2"; shift 2 ;;
      --min-mapq)       MIN_MAPQ="$2"; shift 2 ;;
      --keep-mt)        DROP_MT=0; shift ;;
      --nfr-max)        NFR_MAX="$2"; shift 2 ;;
      --mono-max)       MONO_MAX="$2"; shift 2 ;;
      --macs2-qvalue)   MACS2_QVALUE="$2"; shift 2 ;;
      --toppeaks)       TOPPEAKS_N="$2"; shift 2 ;;
      --igv-region)     IGV_REGION="$2"; shift 2 ;;
      --stages)         STAGES="$2"; shift 2 ;;
      --list-stages)    LIST_STAGES=1; shift ;;
      --force)          FORCE=1; shift ;;
      --check)          CHECK_ONLY=1; shift ;;
      --dry-run)        DRY_RUN=1; shift ;;
      --keep-tmp)       KEEP_TMP=1; shift ;;
      --no-motif)       RUN_MOTIF=0; shift ;;
      --no-atacqc)      RUN_ATACQC=0; shift ;;
      --meme)           RUN_MEME=1; shift ;;
      --tobias)         RUN_TOBIAS=1; shift ;;
      -C|--config)      CONFIG_FILE="$2"; shift 2 ;;
      -h|--help)        usage; exit 0 ;;
      --version)        printf '%s %s\n' "$PROG" "$VERSION"; exit 0 ;;
      *)                usage >&2; die "未知参数: $1" ;;
    esac
  done
}
