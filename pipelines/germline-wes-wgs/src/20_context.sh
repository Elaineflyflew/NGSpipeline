#!/usr/bin/env bash
# =============================================================================
#  src/20_context.sh —— 上下文初始化
#
#  职责：参数校验 → 派生样本名 → 展开 caller/annotator → 建目录 → 重定向日志
#  这是唯一会"改全局状态"的模块，后面的 stage 只读不写。
# =============================================================================

init_context() {
  # ---- 0. 额外配置文件（命令行最后说了算，所以这里只补空缺）----
  if [[ -n "$CONFIG_FILE" ]]; then
    [[ -s "$CONFIG_FILE" ]] || die "配置文件不存在: $CONFIG_FILE"
    # shellcheck source=/dev/null
    source "$CONFIG_FILE"
    log "已加载配置文件: $CONFIG_FILE"
  fi

  # ---- 1. 输入校验 ----
  if [[ -z "$BAM_IN" && ( -z "$R1" || -z "$R2" ) ]]; then
    usage >&2
    die "必须提供 --r1/--r2 或 --bam"
  fi

  if [[ -n "$BAM_IN" ]]; then
    (( DRY_RUN )) || [[ -s "$BAM_IN" ]] || die "BAM 不存在: $BAM_IN"
    [[ -n "$SAMPLE" ]] || SAMPLE=$(basename "$BAM_IN" | sed -E 's/\.(mdup\.|markdup\.|sorted\.|rg\.)?bam$//')
  else
    (( DRY_RUN )) || [[ -s "$R1" ]] || die "R1 不存在: $R1"
    (( DRY_RUN )) || [[ -s "$R2" ]] || die "R2 不存在: $R2"
    [[ -n "$SAMPLE" ]] || SAMPLE=$(basename "$R1" | sed -E 's/(_R?1|_1)?\.(fastq|fq)(\.gz)?$//')
  fi
  [[ -n "$SAMPLE" ]] || die "无法从文件名推导样本名，请用 -s 指定"

  # 参考基因组：dry-run 下允许不存在（方便在笔记本上预览命令）
  if [[ -z "$REF" ]]; then
    die "找不到参考基因组。先跑 deploy/deploy.sh ref，或用 -r 指定"
  fi
  if (( ! DRY_RUN )) && [[ ! -s "$REF" ]]; then
    die "参考基因组不存在: ${REF}（先跑 deploy/deploy.sh ref）"
  fi
  [[ -z "$BED" || -s "$BED" ]] || die "BED 不存在: $BED"
  [[ -z "$BED" ]] && warn "未指定 --bed，将做全基因组检测（WES 会慢很多）"

  # ---- 2. caller / annotator 展开成数组 ----
  case "$CALLER" in
    bcftools)    CALLERS=(bcftools) ;;
    gatk)        CALLERS=(gatk) ;;
    deepvariant) CALLERS=(deepvariant) ;;
    all)         CALLERS=(bcftools gatk deepvariant) ;;
    *) die "--caller 只能是 bcftools|gatk|deepvariant|all，收到: $CALLER" ;;
  esac

  case "$ANNOTATOR" in
    vep)        ANNOTATORS=(vep) ;;
    annovar)    ANNOTATORS=(annovar) ;;
    funcotator) ANNOTATORS=(funcotator) ;;
    all)        ANNOTATORS=(vep funcotator annovar) ;;
    none)       ANNOTATORS=() ;;
    *) die "--annotator 只能是 vep|annovar|funcotator|none|all，收到: $ANNOTATOR" ;;
  esac

  # 空数组在 set -u 下展开会报错，统一转成可打印的字符串
  ANNOT_DESC="none"
  (( ${#ANNOTATORS[@]} )) && ANNOT_DESC="${ANNOTATORS[*]}"

  # ---- 3. 默认阶段：有 BAM 就不用跑前三步 ----
  local default_stages
  if [[ -n "$BAM_IN" ]]; then
    default_stages="bam,call,annotate,report"
  else
    default_stages="qc,trim,align,bam,call,annotate,report"
  fi
  STAGES="${STAGES:-$default_stages}"
  # 非法阶段名直接拦掉，免得写错一个字母就静默跳过整段流程
  local s
  for s in ${STAGES//,/ }; do
    case "$s" in
      qc|trim|align|bam|call|annotate|report) ;;
      *) die "未知阶段: ${s}（可用: qc,trim,align,bam,call,annotate,report）" ;;
    esac
  done

  # ---- 4. 目录骨架（沿用课程里的编号习惯）----
  QC="$OUTDIR/01_qc"; CLEAN="$OUTDIR/02_clean"; MAP="$OUTDIR/03_map"
  STAT="$OUTDIR/04_bamstats"; VAR="$OUTDIR/05_variant"; ANN="$OUTDIR/06_annotate"
  REP="$OUTDIR/07_report"; LOGS="$OUTDIR/00_logs"; STATUS="$OUTDIR/.status"; TMPD="$OUTDIR/.tmp"
  mkdir -p "$QC" "$CLEAN" "$MAP" "$STAT" "$VAR" "$ANN" "$REP" "$LOGS" "$STATUS" "$TMPD"
  local c
  for c in "${CALLERS[@]}"; do mkdir -p "$VAR/$c" "$VAR/$c/plot"; done

  # ---- 5. 派生路径 ----
  MDUP="$MAP/${SAMPLE}.mdup.bam"
  MDUP_BAI="${MDUP}.bai"

  # ---- 6. PATH：项目自带工具优先，其次课程服务器的工具目录 ----
  local p
  for p in "$PROJECT_ROOT/extratools/bin" ${EXTRA_PATH//:/ }; do
    [[ -d "$p" ]] && [[ ":$PATH:" != *":$p:"* ]] && export PATH="$p:$PATH"
  done

  # ---- 7. 日志：dry-run / --check 时留在屏幕，否则落文件 ----
  LOG="$LOGS/${SAMPLE}.$(date +%Y%m%d_%H%M%S).log"
  if (( DRY_RUN || CHECK_ONLY )); then
    :
  elif [[ -t 1 ]]; then
    # 终端里：同时打屏和落盘
    exec > >(tee -a "$LOG") 2>&1
  else
    # nohup / cron：直接写文件。这里不能用 tee，SSH 断开会 SIGPIPE 把主进程带崩
    exec >> "$LOG" 2>&1
  fi
}

# 启动横幅
header() {
  printf '\n=== %s %s | 样本 %s ===\n' "$PROG" "$VERSION" "$SAMPLE"
  printf '参考: %s\n输出: %s\n日志: %s\n\n' "$REF" "$OUTDIR" "$LOG"
}

# 每个 caller 最终 VCF 的路径；注释和报告都走这个函数，避免路径散落
final_vcf() {
  case "$1" in
    bcftools)    printf '%s' "$VAR/bcftools/${SAMPLE}.bcftools.final.vcf.gz" ;;
    gatk)        if [[ -s "$VAR/gatk/${SAMPLE}.gatk.final.vcf.gz" ]]; then
                   printf '%s' "$VAR/gatk/${SAMPLE}.gatk.final.vcf.gz"
                 else
                   printf '%s' "$VAR/gatk/${SAMPLE}.HC.cnn.vcf.gz"
                 fi ;;
    deepvariant) printf '%s' "$VAR/deepvariant/${SAMPLE}.dv.vcf.gz" ;;
  esac
}
