#!/usr/bin/env bash
# =============================================================================
#  src/10_cli.sh —— 命令行接口
#
#  只做两件事：打印帮助（usage）、把 argv 解析成全局变量（parse_cli）。
#  新增一个参数 = 在 usage 里加一行 + 在 parse_cli 的 case 里加一个分支。
# =============================================================================

# 全局状态（由 parse_cli 填充）
R1=""; R2=""; BAM_IN=""; SAMPLE=""
STAGES=""
FORCE=0; DRY_RUN=0; CHECK_ONLY=0; KEEP_TMP=0; LIST_STAGES=0
CONFIG_FILE=""

usage() {
cat <<EOF
${PROG} ${VERSION} —— WES/WGS germline 变异检测流水线

用法
  ${PROG} -1 R1.fq.gz -2 R2.fq.gz [选项]        # 从原始 FASTQ 开始
  ${PROG} --bam sample.mdup.bam   [选项]        # 复用已有 BAM

输入（二选一）
  -1, --r1 FILE          双端 R1（.fastq/.fq，可 .gz）
  -2, --r2 FILE          双端 R2
      --bam FILE         已有 BAM，跳过质控/修剪/比对

核心参数
  -r, --ref FILE         参考基因组 fasta（默认: reference/genome/ 下自动探测）
  -o, --outdir DIR       输出目录（默认: ./results）
  -s, --sample NAME      样本名（默认从文件名推导）
  -t, --threads N        线程数（默认: ${THREADS}）
  -b, --bed FILE         WES 捕获区间 BED（强烈建议，能大幅提速）
  -c, --caller X         bcftools | gatk | deepvariant | all（默认: ${CALLER}）
  -a, --annotator X      vep | annovar | funcotator | none | all（默认: ${ANNOTATOR}）

数据库（默认都指向项目内目录，deploy 之后不用改）
      --bundle DIR           GATK bundle（dbsnp / hapmap / Mills / 1000G）
      --vep-cache DIR        VEP cache（父目录或 homo_sapiens/<ver> 都行）
      --vep-assembly STR     默认 ${VEP_ASSEMBLY}
      --annovar-db DIR       ANNOVAR humandb
      --annovar-protocols S  默认 ${ANNOVAR_PROTOCOLS}
      --annovar-ops S        默认 ${ANNOVAR_OPS}
      --funcotator-data DIR  Funcotator 数据源目录
      --dv-sif FILE          DeepVariant singularity 镜像
      --dv-bind DIR          singularity 绑定路径（默认项目根目录）
      --dv-model X           WES | WGS（默认 ${DV_MODEL}）

阈值
      --min-mapq N           比对质量下限（默认 ${MIN_MAPQ}）
      --min-qual N           bcftools QUAL 下限（默认 ${MIN_QUAL}）
      --min-dp N             bcftools DP 下限（默认 ${MIN_DP}）
      --java-mem X           GATK 堆内存（默认 ${JAVA_MEM}）

流程控制
      --stages LIST          qc,trim,align,bam,call,annotate,report（默认全跑）
      --list-stages          打印所有阶段及说明后退出
      --force                忽略完成标记，强制重跑
      --check                只检查工具和数据库，不跑数据
      --dry-run              只打印命令，不执行
      --keep-tmp             保留中间 BAM
  -C, --config FILE          额外加载一个配置文件（覆盖 config/pipeline.conf）
  -h, --help                 显示本帮助
      --version              显示版本

IGV
      --igv-region STR       截图区域（默认 ${IGV_REGION}，即 TP53）
      --igv-max N            候选位点最多截图数（默认 ${IGV_MAX}）

环境变量
  上面每个参数都有同名大写变量（REF / THREADS / CALLER / BUNDLE ...），
  可以写进 ~/.bashrc 或 env 文件。

示例
  # 全流程 + bcftools + VEP
  ${PROG} -1 raw/Case1_R1.fastq.gz -2 raw/Case1_R2.fastq.gz -c bcftools -a vep -t 16

  # 复用 BAM，GATK + Funcotator
  ${PROG} --bam out/03_map/Case1.mdup.bam -c gatk -a funcotator

  # 三条路线一起跑，自动出一致性对比
  ${PROG} --bam Case1.mdup.bam -c all -a none

  # 上服务器前先看会跑什么
  ${PROG} ... --dry-run
EOF
}

list_stages() {
cat <<'EOF'
阶段        内容                                              主要产物
qc          FastQC 原始数据质控                              01_qc/*.html
trim        fastp 接头/低质量修剪                            02_clean/*.clean.fq.gz
align       bwa mem + samtools view -q                      03_map/*.q20.bam
bam         fixmate → sort → markdup → flagstat/bedcov       03_map/*.mdup.bam, 04_bamstats/
call        bcftools | gatk | deepvariant（可多选）           05_variant/<caller>/
annotate    vep | annovar | funcotator（可多选）              06_annotate/
report      汇总指标 + 一致性对比 + IGV 批处理脚本            07_report/
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
      -1|--r1)             R1="$2"; shift 2 ;;
      -2|--r2)             R2="$2"; shift 2 ;;
      --bam)               BAM_IN="$2"; shift 2 ;;
      -r|--ref)            REF="$2"; shift 2 ;;
      -o|--outdir)         OUTDIR="$2"; shift 2 ;;
      -s|--sample)         SAMPLE="$2"; shift 2 ;;
      -t|--threads)        THREADS="$2"; shift 2 ;;
      -b|--bed)            BED="$2"; shift 2 ;;
      -c|--caller)         CALLER="$2"; shift 2 ;;
      -a|--annotator)      ANNOTATOR="$2"; shift 2 ;;
      --bundle)            BUNDLE="$2"; shift 2 ;;
      --vep-cache)         VEP_CACHE="$2"; shift 2 ;;
      --vep-assembly)      VEP_ASSEMBLY="$2"; shift 2 ;;
      --annovar-db)        ANNOVAR_DB="$2"; shift 2 ;;
      --annovar-protocols) ANNOVAR_PROTOCOLS="$2"; shift 2 ;;
      --annovar-ops)       ANNOVAR_OPS="$2"; shift 2 ;;
      --funcotator-data)   FUNCOTATOR_DATA="$2"; shift 2 ;;
      --dv-sif)            DV_SIF="$2"; shift 2 ;;
      --dv-bind)           DV_BIND="$2"; shift 2 ;;
      --dv-model)          DV_MODEL="$2"; shift 2 ;;
      --min-mapq)          MIN_MAPQ="$2"; shift 2 ;;
      --min-qual)          MIN_QUAL="$2"; shift 2 ;;
      --min-dp)            MIN_DP="$2"; shift 2 ;;
      --java-mem)          JAVA_MEM="$2"; shift 2 ;;
      --stages)            STAGES="$2"; shift 2 ;;
      --igv-region)        IGV_REGION="$2"; shift 2 ;;
      --igv-max)           IGV_MAX="$2"; shift 2 ;;
      -C|--config)         CONFIG_FILE="$2"; shift 2 ;;
      --force)             FORCE=1; shift ;;
      --check)             CHECK_ONLY=1; shift ;;
      --dry-run)           DRY_RUN=1; shift ;;
      --keep-tmp)          KEEP_TMP=1; shift ;;
      --list-stages)       LIST_STAGES=1; shift ;;
      -h|--help)           usage; exit 0 ;;
      --version)           printf '%s %s\n' "$PROG" "$VERSION"; exit 0 ;;
      *) printf '未知参数: %s\n\n' "$1" >&2; usage >&2; exit 2 ;;
    esac
  done
}
