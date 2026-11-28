#!/usr/bin/env bash
# =============================================================================
#  src/31_stage_trim.sh —— ② 接头修剪 + 质量过滤
#
#  这一步在 ATAC-seq 里是可选的：SRA 下来的数据多数已经清洗过，
#  接头残留很少。但保留它有两个好处：
#    · 万一数据没洗过，接头会显著拉低比对率
#    · fastp 的 duplication / adapter 统计本身是一份有用的 QC
#
#  不想跑就在 --stages 里去掉 trim；③ 会自动回退到原始 FASTQ。
# =============================================================================

stage_trim() {
  section "② 修剪 fastp"
  need fastp

  local c1="$CLEAN/${SAMPLE}_R1.clean.fq.gz"
  local c2="$CLEAN/${SAMPLE}_R2.clean.fq.gz"

  run fastp -i "$R1" -I "$R2" -o "$c1" -O "$c2" \
      --detect_adapter_for_pe \
      --qualified_quality_phred 15 \
      --length_required 36 \
      -j "$CLEAN/${SAMPLE}.fastp.json" \
      -h "$CLEAN/${SAMPLE}.fastp.html" \
      --thread "$THREADS"

  if (( ! DRY_RUN )); then
    local before after
    before=$(zcat "$R1" | awk 'END{print NR/4}')
    after=$(zcat "$c1"  | awk 'END{print NR/4}')
    {
      printf 'reads_before\t%s\n' "$before"
      printf 'reads_after\t%s\n'  "$after"
      printf 'retention_pct\t%s\n' "$(pct "$after" "$before")"
    } > "$CLEAN/${SAMPLE}.retention.txt"
    log "保留率: $(grep retention "$CLEAN/${SAMPLE}.retention.txt" | cut -f2)%"
  fi

  mark_done trim
}
