#!/usr/bin/env bash
# =============================================================================
#  src/51_stage_bigwig.sh —— ⑪ bigWig 轨道 + TSS 富集曲线
#
#  bamCoverage 的几个参数不是可选的：
#   --ignoreDuplicates      去掉 PCR 重复，否则高深度区会被重复 read 虚高
#   --blackListFileName     去掉 ENCODE 黑名单（着丝粒、rRNA 等假信号区）
#   --normalizeUsing CPM    不同样本的测序深度不同，不归一化就没法横向比
#   --extendReads           PE 数据按实际片段长度延伸，单端才用 --extendReads N
#   -bs 1                   1 bp 分辨率，IGV 里放大看足迹需要
#
#  TSS 富集曲线是判断"信号是否真的集中在启动子"的直观证据，
#  和 ⑩ 算出来的 TSSE 数值应该互相印证：数值 ≥5 的样本，曲线上 TSS 处应有明显尖峰。
# =============================================================================

stage_bigwig() {
  section "⑪ 生成 bigWig"
  need bamCoverage

  checkfile "$SRT2" "⑥ 产出的 srt2.bam"

  local blackopt=()
  [[ -n "$BLACKLIST" && -s "$BLACKLIST" ]] && blackopt=(--blackListFileName "$BLACKLIST")

  run bamCoverage -p "$THREADS" --ignoreDuplicates \
      ${blackopt[@]+"${blackopt[@]}"} \
      -bs "$BIGWIG_BINSIZE" \
      -b "$SRT2" \
      -o "$BW" \
      --normalizeUsing "$BIGWIG_NORMALIZE" \
      --extendReads

  # ---- TSS 富集曲线 ----
  if have computeMatrix; then
    local tssbed="$DEEP/${SAMPLE}.tss.bed"
    if (( ! DRY_RUN )); then
      # reference-point 只要 TSS 单点，取每个基因的 TSS 位置
      awk -F'\t' '$3=="gene" {
          if (match($0, /gene_type "[^"]+"/)) {
            gt = substr($0, RSTART + 11, RLENGTH - 12)
            if (gt != "protein_coding") next
          }
          tss = ($7 == "+") ? $4 : $5
          if (tss - 1 < 0) next
          print $1 "\t" (tss - 1) "\t" tss "\t" $7
        }' "$GTF" | sort -k1,1 -k2,2n > "$tssbed"
    fi

    run computeMatrix reference-point \
        -S "$BW" -R "$tssbed" \
        -a "$PROFILE_FLANK" -b "$PROFILE_FLANK" \
        -bs "$PROFILE_BINSIZE" -p "$THREADS" \
        -o "$DEEP/${SAMPLE}.tss_matrix.gz" \
        --referencePoint TSS \
        --missingDataAsZero \
        --samplesLabel "$SAMPLE"

    run plotProfile -m "$DEEP/${SAMPLE}.tss_matrix.gz" \
        -o "$DEEP/${SAMPLE}.tss_profile.pdf" \
        --plotTitle "TSS enrichment — $SAMPLE" \
        --perGroup

    run plotHeatmap -m "$DEEP/${SAMPLE}.tss_matrix.gz" \
        -o "$DEEP/${SAMPLE}.tss_heatmap.pdf" \
        --colorMap Blues \
        --heatmapHeight 6 --heatmapWidth 4 \
        --regionsLabel "TSS"
  else
    warn "没有 computeMatrix，跳过 TSS 曲线"
  fi

  # ---- plotFingerprint：判断"信号有没有富集"最直观的一张图 ----
  # 曲线明显偏离对角线 = 有富集；贴着对角线 = 基本是背景
  if have plotFingerprint; then
    run plotFingerprint -b "$SRT2" \
        -plot "$DEEP/${SAMPLE}.fingerprint.pdf" \
        --outRawCounts "$DEEP/${SAMPLE}.fingerprint.tab" \
        --labels "$SAMPLE" \
        -p "$THREADS"
  fi

  mark_done bigwig
}
