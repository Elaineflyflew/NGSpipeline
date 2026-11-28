#!/usr/bin/env bash
# =============================================================================
#  src/42_stage_coverage.sh —— ⑪ bigWig 信号轨 + TSS 富集
#
#  bamCoverage 的三个开关都有理由：
#    --ignoreDuplicates        PCR 重复会在信号图上造成"假高峰"
#    --blackListFileName      着丝粒/端粒/rDNA 高拷贝区同理
#    --extendReads            双端数据按真实片段延伸，这才是 ATAC 的信号单位
#
#  TSS 富集分数 = TSS ±50 bp 均值 / 两端 flank 均值，和 ATACseqQC 的
#  TSSEscore 是同一个定义，⑭ 可以拿来做交叉验证。
# =============================================================================

stage_coverage() {
  section "⑪ bigWig + TSS 富集"
  need bamCoverage
  need computeMatrix
  need plotProfile
  need_file "$BAM_SRT2" "⑥ 的 srt2.bam（先跑 --stages bam）"

  local bw="$COV/${SAMPLE}.bamCoverage.bw"
  local bl=() ig=()
  [[ -n "$BLACKLIST" ]] && bl=(--blackListFileName "$BLACKLIST")
  [[ -n "$IGNORE_FOR_NORM" ]] && ig=(--ignoreForNormalization "$IGNORE_FOR_NORM")

  run bamCoverage -p "$THREADS" --ignoreDuplicates \
      ${bl[@]+"${bl[@]}"} ${ig[@]+"${ig[@]}"} \
      -bs "$BW_BINSIZE" --normalizeUsing "$BW_NORM" --extendReads \
      -b "$BAM_SRT2" -o "$bw"

  # ---- TSS BED：全样本共用，做一次就够 ----
  if tss_bed_ready; then
    log "TSS BED 已存在，复用: $TSS_BED"
  elif [[ -n "$GTF" ]]; then
    run "$PYTHON" "$PROJECT_ROOT/src/py/tss_prep.py" --gtf "$GTF" --out "$TSS_BED"
  else
    warn "没有 GTF，跳过 TSS 富集（deploy/deploy.sh annotation 可补）"
    mark_done coverage
    return 0
  fi

  local matrix="$TSS/${SAMPLE}.tss_matrix.gz"
  run computeMatrix reference-point -S "$bw" -R "$TSS_BED" \
      -b "$TSS_FLANK" -a "$TSS_FLANK" -bs "$TSS_BIN" \
      --referencePoint TSS --missingDataAsZero \
      -p "$THREADS" -o "$matrix"

  run plotProfile -m "$matrix" -o "$TSS/${SAMPLE}.tss_profile.pdf" \
      --outFileNameData "$TSS/${SAMPLE}.tss_profile_data.tsv" \
      --plotTitle "TSS enrichment ($SAMPLE)" \
      --yAxisLabel "mean coverage (CPM)"

  run "$PYTHON" "$PROJECT_ROOT/src/py/tss_score.py" \
      --matrix "$matrix" --flank "$TSS_FLANK" --window "$TSS_WINDOW" \
      --out "$TSS/${SAMPLE}.tss_score.tsv"

  (( DRY_RUN )) && { mark_done coverage; return 0; }

  local key val
  while IFS=$'\t' read -r key val; do
    case "$key" in
      tss_enrichment_score) metric "11_coverage" "tss_enrichment_score" "$val" "fold" ">6" ;;
      tss_signal_mean)      metric "11_coverage" "tss_signal_mean"      "$val" "" "" ;;
      flank_signal_mean)    metric "11_coverage" "flank_signal_mean"    "$val" "" "" ;;
      n_tss_regions)        metric "11_coverage" "n_tss_regions"        "$val" "genes" "" ;;
      *) ;;
    esac
  done < "$TSS/${SAMPLE}.tss_score.tsv"

  local s
  s=$(awk -F'\t' '$1=="tss_enrichment_score"{print $2}' "$TSS/${SAMPLE}.tss_score.tsv")
  ok "TSS 富集倍数 ${s}×（ENCODE 标准 >6）"
  mark_done coverage
}
