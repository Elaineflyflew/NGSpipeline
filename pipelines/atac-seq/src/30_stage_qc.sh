#!/usr/bin/env bash
# =============================================================================
#  src/30_stage_qc.sh —— ① 原始数据质控
#
#  看什么：reads 数、GC 分布、per-base quality、adapter content
#  结论用来定 ② 的参数：有接头才开 --detect_adapter_for_pe，Q30 低就压 -q
# =============================================================================

stage_qc() {
  section "① FastQC 原始质检"
  need fastqc
  fastqc_probe

  run fastqc -t "$THREADS" ${FASTQC_EXTRA[@]+"${FASTQC_EXTRA[@]}"} --outdir "$QC" "$R1" "$R2"
  (( DRY_RUN )) && { mark_done qc; return 0; }

  local stem z n gc
  stem=$(fastqc_stem "$R1")
  z="$QC/${stem}_fastqc.zip"
  n=$(fastqc_field "$z" "Total Sequences" || true)
  gc=$(fastqc_field "$z" "total %GC" || true)

  [[ -n "$n" ]]  && metric "01_fastqc" "total_sequences" "$n"  "reads" ""
  [[ -n "$gc" ]] && metric "01_fastqc" "gc_content"      "$gc" "%"     "38-42"

  ok "FastQC 报告: $QC/*_fastqc.html（reads ${n:-?}，GC ${gc:-?}%）"
  mark_done qc
}
