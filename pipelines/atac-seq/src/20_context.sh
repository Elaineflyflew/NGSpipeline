#!/usr/bin/env bash
# =============================================================================
#  src/20_context.sh —— 上下文初始化
#
#  职责：参数校验 → 派生样本名 → 建目录 → 派生路径 → 重定向日志
#  这是唯一会"改全局状态"的模块，后面的 stage 只读不写。
# =============================================================================

ALL_STAGES="qc,trim,align,bam,shift,fraction,peak,metrics,bigwig,annotate,atacqc,report"

init_context() {
  # ---- 0. 额外配置文件（命令行最后说了算，这里只补空缺）----
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
    [[ -n "$SAMPLE_OPT" ]] || SAMPLE_OPT=$(basename "$BAM_IN" | sed -E 's/\.(srt2|srt|q20|shift)?\.?bam$//')
  else
    (( DRY_RUN )) || [[ -s "$R1" ]] || die "R1 不存在: $R1"
    (( DRY_RUN )) || [[ -s "$R2" ]] || die "R2 不存在: $R2"
    [[ -n "$SAMPLE_OPT" ]] || SAMPLE_OPT=$(basename "$R1" | sed -E 's/(_R?1|_1)?\.(fastq|fq)(\.gz)?$//')
  fi
  SAMPLE="$SAMPLE_OPT"
  [[ -n "$SAMPLE" ]] || die "无法从文件名推导样本名，请用 -s 指定"

  # ---- 2. 参考数据 ----
  (( DRY_RUN )) || checkfile "$BT2_INDEX.1.bt2" "bowtie2 索引"
  (( DRY_RUN )) || checkfile "$BLACKLIST" "ENCODE blacklist"
  (( DRY_RUN )) || checkfile "$GTF" "基因注释 GTF"
  [[ -n "$FASTA" ]] || warn "没找到参考 fasta，bedtools getfasta / TOBIAS 相关步骤会跳过"

  # ---- 3. 阶段白名单 ----
  local default_stages s
  if [[ -n "$BAM_IN" ]]; then
    default_stages="peak,metrics,bigwig,annotate,atacqc,report"
  else
    default_stages="$ALL_STAGES"
  fi
  STAGES="${STAGES:-$default_stages}"
  for s in ${STAGES//,/ }; do
    case ",$ALL_STAGES," in
      *",$s,"*) ;;
      *) die "未知阶段: ${s}（--list-stages 看全部）" ;;
    esac
  done

  # ---- 4. 目录骨架（编号与流程图的 15 步对应）----
  QC="$OUTDIR/01_qc";          CLEAN="$OUTDIR/02_clean";   MAP="$OUTDIR/03_map"
  BAM="$OUTDIR/04_bamfiles";   PEAK="$OUTDIR/05_peaks";    DEEP="$OUTDIR/06_deeptools"
  ANN="$OUTDIR/07_annotation"; MOTIF="$OUTDIR/08_motif";   AQC="$OUTDIR/09_atacqc"
  REP="$OUTDIR/10_report";     LOGS="$OUTDIR/00_logs"
  STATUS="$OUTDIR/.status";    TMPD="$OUTDIR/.tmp"
  mkdir -p "$QC" "$CLEAN" "$MAP" "$BAM" "$PEAK" "$DEEP" "$ANN" "$MOTIF" \
           "$AQC" "$REP" "$LOGS" "$STATUS" "$TMPD"

  # ---- 5. 派生路径（后面所有模块都走这几个变量，不各写各的）----
  SAM="$MAP/${SAMPLE}_bowtie2.sam"
  Q20="$BAM/${SAMPLE}.q20.bam"
  SRT="$BAM/${SAMPLE}.srt.bam"          # 未偏移，给 ATACseqQC / TOBIAS ATACorrect
  SHIFT_BAM="$BAM/${SAMPLE}.shift.bam"
  SRT2="$BAM/${SAMPLE}.srt2.bam"        # 已偏移，peak calling / FRiP / bigWig 都用它
  NFR="$BAM/${SAMPLE}.NFR.bam"
  MONO="$BAM/${SAMPLE}.mono.bam"
  MULTI="$BAM/${SAMPLE}.multi.bam"
  BW="$OUTDIR/07_bigwigs/${SAMPLE}.bamCoverage.bw"
  PEAKS_RAW="$PEAK/${SAMPLE}_peaks.narrowPeak"
  PEAKS="$PEAK/${SAMPLE}.filtered.peaks.narrowPeak"
  PEAKS_NFR="$PEAK/${SAMPLE}_NFR.filtered.peaks.narrowPeak"
  PEAKS_MONO="$PEAK/${SAMPLE}_mono.filtered.peaks.narrowPeak"
  TOPPEAKS="$PEAK/${SAMPLE}.toppeaks.narrowPeak"

  mkdir -p "$OUTDIR/07_bigwigs"

  # --bam 传进来的文件直接接管对应路径，避免下游去找一个不存在的 results 文件
  if [[ -n "$BAM_IN" ]]; then
    if is_shifted "$BAM_IN"; then
      SRT2="$BAM_IN"
    else
      SRT="$BAM_IN"
      warn "输入的 BAM 看起来未做 Tn5 偏移；如需 ⑤，把 shift 加进 --stages"
    fi
  fi

  # ---- 6. PATH：项目自带工具优先，其次课程服务器的工具目录 ----
  local p
  for p in "$PROJECT_ROOT/extratools/bin" ${EXTRA_PATH//:/ }; do
    [[ -d "$p" ]] && [[ ":$PATH:" != *":$p:"* ]] && export PATH="$p:$PATH"
  done

  # ---- 7. 日志 ----
  LOG="$LOGS/${SAMPLE}.$(date +%Y%m%d_%H%M%S).log"
  if (( DRY_RUN || CHECK_ONLY )); then
    :
  elif [[ -t 1 ]]; then
    exec > >(tee -a "$LOG") 2>&1
  else
    # nohup / cron：直接写文件。这里不能用 tee，SSH 断开会 SIGPIPE 把主进程带崩
    exec >> "$LOG" 2>&1
  fi
}

header() {
  printf '\n=== %s %s | 样本 %s ===\n' "$PROG" "$VERSION" "$SAMPLE"
  printf '物种: %s\n索引: %s\n输出: %s\n日志: %s\n\n' \
    "$GENOME" "$BT2_INDEX" "$OUTDIR" "$LOG"
}

# ---------------------------------------------------------------- 路径查询
# 注释和报告都走这几个函数取路径，避免路径散落在各个 stage 里
bam_for_stage() {
  case "$1" in
    srt)  printf '%s' "$SRT" ;;
    srt2) printf '%s' "$SRT2" ;;
    *)    printf '%s' "$SRT2" ;;
  esac
}

peak_for() {
  case "$1" in
    all)  printf '%s' "$PEAKS" ;;
    nfr)  printf '%s' "$PEAKS_NFR" ;;
    mono) printf '%s' "$PEAKS_MONO" ;;
  esac
}

# 是否已经偏移过（--bam 传进来的可能已经是 srt2）
is_shifted() {
  [[ "$1" == *srt2* || "$1" == *shift* ]]
}
