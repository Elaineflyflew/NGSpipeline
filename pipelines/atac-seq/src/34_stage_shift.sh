#!/usr/bin/env bash
# =============================================================================
#  src/34_stage_shift.sh —— ⑤⑥ Tn5 偏移校正
#
#  原理：Tn5 转座酶切开 DNA 后会留下 9 bp 的错位缺口（正链 +4，负链 -5）。
#  直接把 read 的 5' 端当成插入位点会整体偏 4-5 bp，在核小体边界和
#  footprint 分析里这个偏差足以让结论反过来。
#
#  alignmentSieve --ATACshift 做的事就是把每条 read 的比对位置按链方向挪一下。
#
#  坑：alignmentSieve 输出的 BAM 只保留前 11 个字段（丢掉后面的 tag），
#  而且不再是坐标有序的，所以必须重新 sort + index 才能给 MACS2 用。
# =============================================================================

stage_shift() {
  section "⑤ Tn5 偏移校正 alignmentSieve"
  need alignmentSieve
  need samtools

  if (( ! ATAC_SHIFT )); then
    warn "ATAC_SHIFT=0，跳过偏移；直接用未偏移 BAM（footprint 结果会有 4-5 bp 偏差）"
    SRT2="$SRT"
    mark_done shift
    return 0
  fi

  if [[ -n "$BAM_IN" ]] && is_shifted "$BAM_IN"; then
    log "输入 BAM 已偏移，跳过 ⑤"
    mark_done shift
    return 0
  fi

  checkfile "$SRT" "④ 产出的 srt.bam"

  run alignmentSieve -p "$THREADS" --ATACshift -b "$SRT" -o "$SHIFT_BAM"

  section "⑥ 重新排序建索引"
  run samtools sort -@ "$THREADS" -m 1G -T "$TMPD/srt2.$SAMPLE" -o "$SRT2" "$SHIFT_BAM"
  run samtools index -@ "$THREADS" "$SRT2"

  if (( ! DRY_RUN )); then
    samtools flagstat -@ "$THREADS" "$SRT2" > "$BAM/${SAMPLE}.srt2.flagstat.txt"
    log "已偏移 BAM: $SRT2"
  fi

  (( KEEP_TMP )) || { (( DRY_RUN )) || rm -f "$SHIFT_BAM"; }
  mark_done shift
}
