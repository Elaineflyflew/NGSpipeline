#!/usr/bin/env bash
# =============================================================================
#  src/20_context.sh —— 上下文初始化
#
#  职责：加载额外配置 → 校验参数 → 派生样本名与路径 → 建目录 → 改 PATH → 落日志
#  这是唯一会"改全局状态"的模块，后面的 stage 只读不写。
# =============================================================================

ALL_STAGES="qc,trim,align,bam,frag,peaks,frip,coverage,annot,motif,atacqc,igv,report"

init_context() {
  # ---- 0. 额外配置文件（命令行最后说了算，所以这里只补空缺）----
  if [[ -n "$CONFIG_FILE" ]]; then
    [[ -s "$CONFIG_FILE" ]] || die "配置文件不存在: $CONFIG_FILE"
    # shellcheck source=/dev/null
    source "$CONFIG_FILE"
    log "已加载配置文件: $CONFIG_FILE"
  fi

  # ---- 1. 输入校验 ----
  # --check 只查工具和参考数据，不碰输入，所以允许不给 -1/-2
  if (( ! CHECK_ONLY )); then
    if [[ -z "$R1" || -z "$R2" ]]; then
      usage >&2
      die "必须提供 -1/--r1 和 -2/--r2"
    fi
    (( DRY_RUN )) || [[ -s "$R1" ]] || die "R1 不存在: $R1"
    (( DRY_RUN )) || [[ -s "$R2" ]] || die "R2 不存在: $R2"
  fi

  # sub20M.SRR13849313_1.fastq.gz → sub20M.SRR13849313
  if [[ -z "$SAMPLE" ]]; then
    SAMPLE=$(basename "${R1:-sample}" | sed -E 's/(_R?1|_1)?\.(fastq|fq)(\.gz)?$//')
    [[ -n "$SAMPLE" ]] || SAMPLE="sample"
  fi

  # ---- 2. 参考数据：空则自动探测，dry-run 下允许不存在 ----
  [[ -n "$GENOME" ]] || GENOME=$(discover_genome)
  [[ -n "$BT2" ]]    || BT2=$(discover_bt2)
  [[ -n "$GTF" ]]    || GTF=$(discover_one "$PROJECT_ROOT/reference/annotation" '*.gtf' '*.gtf.gz')
  [[ -n "$BLACKLIST" ]] || BLACKLIST=$(discover_one "$PROJECT_ROOT/reference/blacklist" '*.bed' '*.bed.gz')

  # 基因组缺不缺，交给 preflight 判断——因为只有 align 之后的阶段才用得上它，
  # 只跑 --stages qc,trim 的时候不该被"没有基因组"拦住。

  # ---- 3. 阶段列表校验：写错一个字母不该静默跳过整段流程 ----
  STAGES="${STAGES:-$ALL_STAGES}"
  local s
  for s in ${STAGES//,/ }; do
    case ",$ALL_STAGES," in
      *",$s,"*) ;;
      *) die "未知阶段: ${s}（可用: ${ALL_STAGES}）" ;;
    esac
  done

  # ---- 4. 目录骨架（编号 = 作业图上的步骤号）----
  QC="$OUTDIR/01_fastqc";  QCC="$QC/clean"
  TRIM="$OUTDIR/02_fastp"; MAP="$OUTDIR/03_map"; BAM="$OUTDIR/04_bam"
  FRAG="$OUTDIR/05_fragbin"; PEAK="$OUTDIR/06_peaks"; FRIP="$OUTDIR/07_frip"
  COV="$OUTDIR/08_bigwig"; TSS="$OUTDIR/09_tss"; ANN="$OUTDIR/10_annot"
  MOTIF="$OUTDIR/11_motif"; ATACQC="$OUTDIR/12_atacqc"; IGV="$OUTDIR/13_igv"
  METRIC="$OUTDIR/metrics"; REP="$OUTDIR/99_report"
  LOGS="$OUTDIR/00_logs"; STATUS="$OUTDIR/.status"; TMPD="$OUTDIR/.tmp"
  mkdir -p "$QC" "$QCC" "$TRIM" "$MAP" "$BAM" "$FRAG" "$PEAK" "$FRIP" \
           "$COV" "$TSS" "$ANN" "$MOTIF" "$ATACQC" "$IGV" "$METRIC" "$REP" \
           "$LOGS" "$STATUS" "$TMPD"

  # ---- 5. 派生路径（全项目只在这里出现一次）----
  BAM_SRT="$BAM/${SAMPLE}.srt.bam"
  BAM_SHIFT="$BAM/${SAMPLE}.shift.bam"
  BAM_SRT2="$BAM/${SAMPLE}.srt2.bam"
  BAM_NFR="$BAM/${SAMPLE}_nucleosome_free.bam"
  BAM_MONO="$BAM/${SAMPLE}_mono_nucleosome.bam"
  BAM_POLY="$BAM/${SAMPLE}_multi_nucleosome.bam"
  TSS_BED="$TSS/genes.tss.bed"
  FRACS="all nfr mono"

  # ---- 6. 文件句柄上限 ----
  # deeptools 的 alignmentSieve / samtools sort 会同时开很多临时文件，
  # 默认 1024 直接报 "Too many open files"（而且报的是 /tmp 里的临时 BAM）。
  # 软上限能提就提，提不动就算了，后面 sieve 会用更少的线程兜住。
  local _cur _hard
  _cur=$(ulimit -n 2>/dev/null || echo 1024)
  if [[ "$_cur" != "unlimited" ]]; then
    # 不要写死 65536：服务器的硬上限常常正好是 65535，多要一个就整个失败，
    # 于是每次运行都刷一条其实无害的 WARN。先问硬上限，再提到硬上限。
    _hard=$(ulimit -Hn 2>/dev/null || echo "$_cur")
    if [[ "$_hard" == "unlimited" ]] || { [[ "$_hard" =~ ^[0-9]+$ ]] && (( _hard > _cur )); }; then
      ulimit -n "$_hard" 2>/dev/null || true
      _cur=$(ulimit -n 2>/dev/null || echo "$_cur")
    fi
    # 真正会炸的量级是 1024（deeptools 每线程一个临时 BAM），4096 以上都够用。
    if [[ "$_cur" != "unlimited" ]] && (( _cur < 4096 )); then
      warn "文件句柄上限只有 ${_cur}，alignmentSieve 可能报 Too many open files（跑前先 ulimit -n 8192）"
    fi
  fi

  # ---- 7. PATH：项目自带工具优先，其次课程服务器的工具目录 ----
  local p
  for p in "$PROJECT_ROOT/extratools/bin" ${EXTRA_PATH//:/ }; do
    [[ -d "$p" ]] && [[ ":$PATH:" != *":$p:"* ]] && export PATH="$p:$PATH"
  done

  # 统一解释器：优先挑带 matplotlib 的（⑦⑫ 要画图）
  PYTHON=$(pick_python)

  # ---- 8. 日志：dry-run / --check 时留在屏幕，否则落文件 ----
  LOG="$LOGS/${SAMPLE}.$(date +%Y%m%d_%H%M%S).log"
  if (( DRY_RUN || CHECK_ONLY )); then
    :
  elif [[ -t 1 ]]; then
    exec > >(tee -a "$LOG") 2>&1      # 终端里：同时打屏和落盘
  else
    exec >> "$LOG" 2>&1               # nohup/cron：直接写文件，不用 tee（SSH 断开会 SIGPIPE）
  fi
}

# 启动横幅
header() {
  printf '\n=== %s %s | 样本 %s ===\n' "$PROG" "$VERSION" "$SAMPLE"
  printf '输入: %s\n      %s\n' "$R1" "$R2"
  printf '参考: %s\n' "${GENOME:-（未找到）}"
  printf '输出: %s\n日志: %s\n\n' "$OUTDIR" "$LOG"
}

# 三条 peak calling 路线各自的 BAM；⑧⑩ 都走这个函数，避免路径散落
frac_bam() {
  case "$1" in
    all)  printf '%s' "$BAM_SRT2" ;;
    nfr)  printf '%s' "$BAM_NFR" ;;
    mono) printf '%s' "$BAM_MONO" ;;
    *)    die "未知路线: $1（可用: all nfr mono）" ;;
  esac
}
