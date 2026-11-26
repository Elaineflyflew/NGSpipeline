#!/usr/bin/env bash
# =============================================================================
#  src/33_stage_bam.sh —— ④ BAM 精修
#
#  顺序不能改：sort -n → fixmate -m → sort（坐标序）→ markdup
#  fixmate 靠 read name 找配对并写 ms/MC tag；顺序错了不报错，
#  但 markdup 的 dup rate 会静默偏低——这是最典型的"看起来跑通了其实错了"。
# =============================================================================

stage_bam() {
  section "④ BAM 精修（fixmate → sort → markdup）"
  need samtools

  local src
  if [[ -n "$BAM_IN" ]]; then src="$BAM_IN"; else src="$MAP/${SAMPLE}.q${MIN_MAPQ}.bam"; fi
  (( DRY_RUN )) || [[ -s "$src" ]] || die "找不到输入 BAM: $src"

  local nsort="$MAP/${SAMPLE}.nsort.bam"
  local fix="$MAP/${SAMPLE}.fixmate.bam"
  local psort="$MAP/${SAMPLE}.psort.bam"

  run samtools sort -n -@ "$THREADS" -m 1G -T "$MAP/tmp.$SAMPLE" -o "$nsort" "$src"
  run samtools fixmate -m -@ "$THREADS" "$nsort" "$fix"
  run samtools sort -@ "$THREADS" -m 1G -T "$MAP/tmp.$SAMPLE" -o "$psort" "$fix"
  # 不加 -r：保留 DUP 标记而不删除。dup rate 本身是质控指标，IGV 里也要能看到
  run samtools markdup -@ "$THREADS" -s -f "$STAT/${SAMPLE}.markdup.txt" "$psort" "$MDUP"
  run samtools index -@ "$THREADS" "$MDUP"

  if (( ! DRY_RUN )); then
    samtools flagstat -@ "$THREADS" "$MDUP" > "$STAT/${SAMPLE}.flagstat.txt"
    samtools stats    -@ "$THREADS" "$MDUP" > "$STAT/${SAMPLE}.stats.txt"
    samtools idxstats "$MDUP"                > "$STAT/${SAMPLE}.idxstats.txt"
    if command -v plot-bamstats >/dev/null 2>&1; then
      mkdir -p "$STAT/plot"
      plot-bamstats -p "$STAT/plot" "$STAT/${SAMPLE}.stats.txt" >/dev/null 2>&1 || true
    fi
    # 有 BED 时算目标区平均深度，比全基因组平均深度有意义得多
    if [[ -n "$BED" ]]; then
      samtools bedcov "$BED" "$MDUP" \
        | awk '{len+=$3-$2; cov+=$5} END{printf "target_bases\t%d\nmean_depth\t%.1f\n", len, (len>0?cov/len:0)}' \
        > "$STAT/${SAMPLE}.target_depth.txt"
      log "目标区深度: $(grep mean_depth "$STAT/${SAMPLE}.target_depth.txt")"
    fi
  fi

  # 中间文件默认清掉：一个 WES 样本的 nsort/fix/psort 加起来约 3 倍 BAM 大小
  if (( ! KEEP_TMP )) && (( ! DRY_RUN )); then rm -f "$nsort" "$fix" "$psort"; fi
  mark_done bam
}
