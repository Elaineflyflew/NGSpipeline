#!/usr/bin/env bash
# =============================================================================
#  src/31_stage_trim.sh —— ② 修剪
#
#  参数来自 lab5 脚本：--detect_adapter_for_pe -l 36
#  -l 36 是下限而不是"剪到 36"：短于 36 的 read 直接丢，比对率反而更高。
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

  # 修剪前后的 read 数对比是判断"接头污染有多严重"的直接证据
  if (( ! DRY_RUN )); then
    local before after
    before=$(zcat "$R1" | awk 'END{print NR/4}')
    after=$(zcat "$c1"  | awk 'END{print NR/4}')
    printf 'reads_before\t%s\nreads_after\t%s\nretention\t%.2f%%\n' \
      "$before" "$after" "$(awk -v a="$after" -v b="$before" 'BEGIN{print (b>0?100*a/b:0)}')" \
      > "$CLEAN/${SAMPLE}.retention.txt"
    log "保留率: $(tail -1 "$CLEAN/${SAMPLE}.retention.txt")"
  fi

  mark_done trim
}
