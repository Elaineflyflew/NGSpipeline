#!/usr/bin/env bash
# =============================================================================
#  src/50_stage_metrics.sh —— ⑩ FRiP 双算法交叉验证 + TSS 富集分数
#
#  FRiP（Fraction of Reads in Peaks）是 ATAC-seq 最核心的质控指标：
#  ENCODE 的底线是 0.2，0.3 以上算好，0.5 以上属于非常干净的数据。
#
#  这里用两种完全独立的算法算同一个数：
#    A. bedtools coverage 逐 peak 数 reads 再求和（纯计数，不依赖 deeptools）
#    B. deeptools plotEnrichment 的 outRawCounts（走另一套计数实现）
#  两者相差 >5% 就说明其中一条路有问题（常见原因：BED 坐标口径不一致、
#  一个用了 -q 一个没用）。这是很划算的交叉验证 —— 多花 1 分钟，
#  能挡住"FRiP 算错但没人发现"这种最尴尬的错误。
#
#  TSS 富集分数（TSSE）：TSS ± 50 bp 的平均深度 / 两侧 ±(1900-2000) 的背景深度。
#  ENCODE 要求 ≥ 5。它比 FRiP 更能反映"信号是否真的富集在调控区"，
#  因为 FRiP 只看 peak 内，而 peak 本身可能是假的。
# =============================================================================

# 从 GTF 生成 TSS 中心窗和背景窗（BED，0-based）
make_tss_bed() {
  local gtf="$1" center="$2" flank="$3"
  awk -F'\t' '
    $3 == "gene" {
      gt = ""
      if (match($0, /gene_type "[^"]+"/)) {
        gt = substr($0, RSTART + 11, RLENGTH - 12)
      }
      if (gt != "" && gt != "protein_coding") next
      tss = ($7 == "+") ? $4 : $5
      if (tss - 2001 < 0) next
      print $1 "\t" (tss - 51)   "\t" (tss + 50)   > center
      print $1 "\t" (tss - 2001) "\t" (tss - 1900) > flank
      print $1 "\t" (tss + 1900) "\t" (tss + 1999) > flank
    }' "$gtf"
}

# 一个 BED 区域上的平均深度
mean_depth() {
  samtools bedcov -Q "$MIN_MAPQ" "$1" "$2" 2>/dev/null \
    | awk '{ len += $3 - $2; cov += $4 } END { printf "%.4f", (len > 0 ? cov / len : 0) }'
}

stage_metrics() {
  section "⑩ FRiP + TSS 富集"
  need bedtools; need samtools

  checkfile "$SRT2"  "⑥ 产出的 srt2.bam"
  checkfile "$PEAKS" "⑨ 过滤后的 peak"

  local M="$AQC"
  mkdir -p "$M"

  # ---- 算法 A：bedtools coverage ----
  runpipe "bedtools coverage -a '$PEAKS' -b '$SRT2' -counts > '$M/${SAMPLE}.peak_coverage.bed'"

  if (( ! DRY_RUN )); then
    local inpeak total
    inpeak=$(awk '{ s += $NF } END { printf "%d", s }' "$M/${SAMPLE}.peak_coverage.bed")
    total=$(samtools view -@ "$THREADS" -c "$SRT2")
    printf 'reads_total\t%s\n'   "$total"  > "$M/${SAMPLE}.frip.tsv"
    printf 'reads_in_peak\t%s\n' "$inpeak" >> "$M/${SAMPLE}.frip.tsv"
    printf 'frip_bedtools\t%s\n' "$(awk -v a="$inpeak" -v b="$total" 'BEGIN{printf "%.4f", (b>0?a/b:0)}')" >> "$M/${SAMPLE}.frip.tsv"
  fi

  # ---- 算法 B：deeptools plotEnrichment ----
  if have plotEnrichment; then
    run plotEnrichment -b "$SRT2" --BED "$PEAKS" \
        --plotFile "$M/${SAMPLE}.FRiP.pdf" \
        --outRawCounts "$M/${SAMPLE}.FRiP_raw.tab" \
        --plotTitle "FRiP" --labels "peaks" -p "$THREADS"
    if (( ! DRY_RUN )); then
      # 表头：file featureType percent featureReadCount totalReadCount
      local p2
      p2=$(awk -F'\t' 'NR==2{print $3}' "$M/${SAMPLE}.FRiP_raw.tab")
      printf 'frip_deeptools\t%s\n' "$(awk -v x="$p2" 'BEGIN{printf "%.4f", x/100}')" \
        >> "$M/${SAMPLE}.frip.tsv"
    fi
  else
    warn "没有 plotEnrichment，跳过 FRiP 算法 B"
  fi

  # ---- TSS 富集 ----
  local center="$M/${SAMPLE}.tss.center.bed" flank="$M/${SAMPLE}.tss.flank.bed"
  if (( ! DRY_RUN )); then
    make_tss_bed "$GTF" "$center" "$flank"
    local cd fd
    cd=$(mean_depth "$center" "$SRT2")
    fd=$(mean_depth "$flank"  "$SRT2")
    {
      printf 'tss_center_depth\t%s\n' "$cd"
      printf 'tss_flank_depth\t%s\n'  "$fd"
      printf 'tss_enrichment\t%s\n' "$(awk -v a="$cd" -v b="$fd" 'BEGIN{printf "%.2f", (b>0?a/b:0)}')"
      printf 'tss_sites\t%s\n' "$(wc -l < "$center" | tr -d ' ')"
    } > "$M/${SAMPLE}.tsse.tsv"
    log "FRiP: $(awk -F'\t' '$1=="frip_bedtools"{print $2}' "$M/${SAMPLE}.frip.tsv") | TSSE: $(awk -F'\t' '$1=="tss_enrichment"{print $2}' "$M/${SAMPLE}.tsse.tsv")"
  fi

  mark_done metrics
}
