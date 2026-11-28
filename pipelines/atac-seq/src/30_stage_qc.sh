#!/usr/bin/env bash
# =============================================================================
#  src/30_stage_qc.sh —— ① 原始数据质控
#
#  看什么：per-base quality、GC 分布、adapter content、duplication
#  决定什么：Q30 占比和接头污染程度决定 ② 的 fastp 参数要不要更激进
#
#  注意：ATAC-seq 的 read 有相当比例落在 chrM 和开放区，GC 分布比 WGS 更宽，
#  不要拿 WGS 的"GC 曲线必须尖"去要求它。
# =============================================================================

# 从 fastqc_data.txt 抽一行摘要。Q30 要从 "Per sequence quality scores" 段算：
# 那一节给的是"每个质量值上有多少条 read"，不是平均质量。
fastqc_row() {
  local zip="$1" base
  base=$(basename "$zip" _fastqc.zip)
  unzip -p "$zip" "${base}_fastqc/fastqc_data.txt" 2>/dev/null \
    | awk -v s="$base" '
        /^Total Sequences/              { tot = $3 }
        /^Sequence length/              { len = $3 }
        /^%GC/                          { gc  = $2 }
        /^>>Per sequence quality scores/ { inq = 1; next }
        inq && /^>>/                    { inq = 0 }
        inq && $1 ~ /^[0-9]+$/          { all += $2; if ($1 + 0 >= 30) q30 += $2 }
        END {
          printf "%s\t%s\t%s\t%s\t%.1f\n", s, tot, len, gc, (all > 0 ? 100 * q30 / all : 0)
        }'
}

stage_qc() {
  section "① 质控 FastQC"
  need fastqc

  run fastqc -t "$THREADS" --outdir "$QC" "$R1" "$R2"

  if (( ! DRY_RUN )); then
    local f
    printf 'file\ttotal_reads\tread_length\tgc_pct\tq30_pct\n' > "$QC/fastqc.summary.tsv"
    for f in "$QC"/*_fastqc.zip; do
      [[ -f "$f" ]] || continue
      fastqc_row "$f" >> "$QC/fastqc.summary.tsv"
    done
    log "FastQC 汇总: $QC/fastqc.summary.tsv"
  fi

  mark_done qc
}
