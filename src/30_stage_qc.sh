#!/usr/bin/env bash
# =============================================================================
#  src/30_stage_qc.sh —— ① 原始数据质控
#
#  指标：per-base quality、GC 分布、adapter content、duplication
#  结论：Q30 占比、是否存在接头污染，决定后面 fastp 的参数
# =============================================================================

stage_qc() {
  section "① 质控 FastQC"
  need fastqc
  run fastqc -t "$THREADS" --outdir "$QC" "$R1" "$R2"
  log "FastQC 报告: $QC/*_fastqc.html"
  mark_done qc
}
