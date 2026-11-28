#!/usr/bin/env bash
# =============================================================================
#  src/41_stage_frip.sh —— ⑩ FRiP + signal-to-noise
#
#  FRiP 算两遍是有意的：
#    A 路  samtools view -c -L        走 htslib 的区间索引
#    B 路  bedtools bamtobed|intersect 走完全不同的区间求交实现
#  两个数一致，说明 BED 与 BAM 的染色体命名、坐标体系、排序状态都对得上。
#  差超过 0.02 就该停下来查文件，而不是解释生物学。
#
#  plotFingerprint 是"peak 是不是真的"最直观的一张图：曲线几乎贴对角线，
#  说明 reads 在基因组上均匀铺开，根本没有富集。
# =============================================================================

stage_frip() {
  section "⑩ FRiP（双算法交叉验证）"
  need samtools
  need bedtools

  local frac bam peaks
  for frac in $FRACS; do
    bam=$(frac_bam "$frac")
    peaks="$PEAK/${SAMPLE}.${frac}.filtered.narrowPeak"
    need_file "$peaks" "⑨ 的 ${frac} peak（先跑 --stages peaks）"
    run "$PYTHON" "$PROJECT_ROOT/src/py/frip.py" \
        --bam "$bam" --peaks "$peaks" --label "$frac" \
        --outdir "$FRIP/$SAMPLE" --bedtools --threads "$THREADS"
  done

  # ---- 把 <frac>.frip.tsv（label/key/value/unit/bench）翻成标准五列指标 ----
  # 主指标统一命名成 frip_<frac>，报告只认这个名字。
  if (( ! DRY_RUN )); then
    for frac in $FRACS; do
      local label key val unit bench
      while IFS=$'\t' read -r label key val unit bench; do
        [[ -z "$key" ]] && continue
        case "$key" in
          frip_samtools) metric "10_frip" "frip_${frac}" "$val" "$unit" "$bench" ;;
          frip_delta)
            if [[ "$frac" == "all" ]]; then
              metric "10_frip" "frip_delta" "$val" "$unit" "$bench"
            else
              metric "10_frip" "frip_delta_${frac}" "$val" "$unit" "$bench"
            fi ;;
          *) metric "10_frip" "${frac}_${key}" "$val" "$unit" "$bench" ;;
        esac
      done < "$FRIP/$SAMPLE/${frac}.frip.tsv"
    done
  fi

  # ---- signal-to-noise ----
  if find_tool plotFingerprint >/dev/null 2>&1 || (( DRY_RUN )); then
    need plotFingerprint
    run plotFingerprint -b "$BAM_SRT2" --labels "$SAMPLE" \
        --plotFile "$FRIP/${SAMPLE}.plotFingerprint.png" \
        --outQualityMetrics "$FRIP/${SAMPLE}.plotFingerprint.qc.tsv" \
        --numberOfSamples 100000 -p "$THREADS"
  else
    warn "没有 plotFingerprint（deeptools），跳过 signal-to-noise 曲线"
  fi

  (( DRY_RUN )) && { mark_done frip; return 0; }

  local f
  f=$(awk -F'\t' '$2=="frip_all"{print $3}' "$METRIC/10_frip.tsv" 2>/dev/null | head -1)
  ok "FRiP（全部片段）= ${f:-?}"
  mark_done frip
}
