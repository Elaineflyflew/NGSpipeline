#!/usr/bin/env bash
# =============================================================================
#  src/40_stage_peaks.sh —— ⑧ MACS2 peak calling  ⑨ 去 ENCODE 黑名单
#
#  三条路线各 call 一次：全部片段 / NFR / 单核小体。
#  多核小体档**不参与**：跨多个核小体的片段两端落在不同开放区，peak 边界
#  会被拉糊。这一档的用途是看核小体周期性（⑦⑭），不是定位开放区。
#
#  -f BAMPE 用 reads 的真实跨度定义片段；此时 MACS2 会把 --shift/--extsize
#  强制归零（OptValidator.py），传了不报错但静默无效，所以不传。
# =============================================================================

stage_peaks() {
  section "⑧ MACS2 peak calling（3 条路线）"
  need macs2

  local frac bam
  for frac in $FRACS; do
    bam=$(frac_bam "$frac")
    need_file "$bam" "⑦ 的 ${frac} 路线 BAM（先跑 --stages frag）"
    log "路线 $frac: $(basename "$bam")"
    run macs2 callpeak -t "$bam" -f BAMPE "${MACS2_OPTS[@]}" \
        -g "$MACS2_GSIZE" -q "$MACS2_QVALUE" \
        -n "${SAMPLE}.${frac}" --outdir "$PEAK" 2> "$PEAK/${SAMPLE}.${frac}.macs2.log"
  done

  # 指标要读 narrowPeak 文件，dry-run 时文件还不存在
  if (( ! DRY_RUN )); then
    for frac in $FRACS; do
      local n
      n=$(wc -l < "$PEAK/${SAMPLE}.${frac}_peaks.narrowPeak" 2>/dev/null || echo 0)
      metric "08_macs2" "peaks_${frac}" "$n" "peaks" "20k-120k"
    done
  fi

  # ---- ⑨ 黑名单 ----
  section "⑨ 去 ENCODE 黑名单"
  need bedtools
  if [[ -z "$BLACKLIST" ]]; then
    warn "没有黑名单 BED，⑨ 退化成原样拷贝（跑 deploy/deploy.sh blacklist 可补）"
  fi

  for frac in $FRACS; do
    local raw="$PEAK/${SAMPLE}.${frac}_peaks.narrowPeak"
    local filt="$PEAK/${SAMPLE}.${frac}.filtered.narrowPeak"
    if [[ -n "$BLACKLIST" ]]; then
      runpipe "bedtools intersect -v -a '$raw' -b '$BLACKLIST' > '$filt'"
    else
      runpipe "cp '$raw' '$filt'"
    fi
  done

  (( DRY_RUN )) && { mark_done peaks; return 0; }

  for frac in $FRACS; do
    local raw="$PEAK/${SAMPLE}.${frac}_peaks.narrowPeak"
    local filt="$PEAK/${SAMPLE}.${frac}.filtered.narrowPeak"
    local before after pct
    before=$(wc -l < "$raw" 2>/dev/null || echo 0)
    after=$(wc -l  < "$filt" 2>/dev/null || echo 0)
    pct=$(awk -v a="$before" -v b="$after" 'BEGIN{printf "%.2f", (a>0? (a-b)*100/a : 0)}')

    metric "09_blacklist" "peaks_${frac}_before_blacklist" "$before" "peaks" ""
    metric "09_blacklist" "peaks_${frac}_after_blacklist"  "$after"  "peaks" "20k-120k"
    metric "09_blacklist" "blacklist_removed_${frac}"      "$pct"    "%"     "<1"
    ok "路线 $frac: $before → $after（剔除 ${pct}%）"
  done

  mark_done peaks
}
