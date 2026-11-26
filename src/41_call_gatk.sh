#!/usr/bin/env bash
# =============================================================================
#  src/41_call_gatk.sh —— GATK 路线（对应 lab7）
#
#  AddOrReplaceReadGroups → BQSR → HaplotypeCaller(-ERC GVCF) → GenotypeGVCFs
#  → NVScoreVariants(CNN) → FilterVariantTranches(99.95/99.4) → Funcotator
#
#  两个"缺数据就降级而不是崩"的地方：
#  · 没有 bundle → 跳过 BQSR（结果可用，只是少了碱基质量校正）
#  · 没有 hapmap → 跳过 Tranche（FILTER 列不分层，但 VCF 仍是合法输出）
# =============================================================================

call_gatk() {
  section "⑤-B GATK 变异检测"
  need gatk; need samtools

  local d="$VAR/gatk"
  local jopt="-Xmx${JAVA_MEM} -Djava.io.tmpdir=$TMPD"
  local bedopt=()
  [[ -n "$BED" ]] && bedopt=(-L "$BED")

  local rg_bam="$d/${SAMPLE}.rg.bam"
  local srt="$d/${SAMPLE}.rg.srt.bam"
  local mdup="$d/${SAMPLE}.mdup.bam"
  local table="$d/${SAMPLE}.recal.table"
  local bqsr="$d/${SAMPLE}.BQSR.bam"
  local gvcf="$d/${SAMPLE}.HC.g.vcf.gz"
  local hc="$d/${SAMPLE}.HC.vcf.gz"
  local cnn="$d/${SAMPLE}.HC.cnn.vcf.gz"
  local fin="$d/${SAMPLE}.gatk.final.vcf.gz"

  # ---- 1. read group：GATK 硬性要求，只在缺失时补 ----
  local hc_in="$MDUP"
  if (( ! DRY_RUN )) && ! bam_has_rg "$MDUP"; then
    log "BAM 缺少 @RG，补 AddOrReplaceReadGroups"
    run gatk AddOrReplaceReadGroups -I "$MDUP" -O "$rg_bam" \
        -ID "${SAMPLE}_Run1" -SM "$SAMPLE" -LB "$SAMPLE" -PL Illumina -PU Run1
    run gatk SortSam -I "$rg_bam" -O "$srt" -SO coordinate --TMP_DIR "$TMPD"
    run gatk MarkDuplicates -I "$srt" -O "$mdup" --TMP_DIR "$TMPD" \
        -M "$d/${SAMPLE}.mdup.metrics.txt" --CREATE_INDEX true
    hc_in="$mdup"
  fi

  # ---- 2. BQSR：需要 dbsnp + Mills + 1000G ----
  local known1 known2 dbsnp
  known1=$(bundle_file 1000G_phase1.snps.high_confidence.hg38 || true)
  known2=$(bundle_file Mills_and_1000G_gold_standard.indels.hg38 || true)
  dbsnp=$(bundle_file dbsnp_146.hg38 || true)

  # 未压缩的 bundle VCF 需要 GATK 的 .idx；缺了就地补建，
  # 否则 BaseRecalibrator 直接报 "An index is required"
  local kf
  for kf in "$known1" "$known2" "$dbsnp"; do
    [[ -n "$kf" ]] && ensure_vcf_index "$kf"
  done

  if [[ -n "$known1" && -n "$known2" && -n "$dbsnp" ]]; then
    run gatk --java-options "$jopt" BaseRecalibrator -R "$REF" -I "$hc_in" \
        --known-sites "$known1" --known-sites "$known2" --known-sites "$dbsnp" \
        -O "$table"
    run gatk --java-options "$jopt" ApplyBQSR -R "$REF" -I "$hc_in" \
        --bqsr-recal-file "$table" -O "$bqsr"
    hc_in="$bqsr"
  else
    warn "GATK bundle 不完整（需要 dbsnp_146 / Mills / 1000G_phase1），跳过 BQSR"
  fi

  # ---- 3. call ----
  run gatk --java-options "$jopt" HaplotypeCaller -R "$REF" -I "$hc_in" -O "$gvcf" \
      -ERC GVCF --native-pair-hmm-threads "$THREADS" ${bedopt[@]+"${bedopt[@]}"}
  run gatk --java-options "$jopt" GenotypeGVCFs -R "$REF" -V "$gvcf" -O "$hc"

  # ---- 4. CNN 打分：需要 GATK 环境里的 python（lab7 的做法）----
  local gatk_py="$PROJECT_ROOT/extratools/bin"
  [[ -d /home/biotools/miniforge3/envs/gatk/bin ]] && gatk_py="/home/biotools/miniforge3/envs/gatk/bin"
  [[ -d "$gatk_py" ]] && [[ ":$PATH:" != *":$gatk_py:"* ]] && export PATH="$gatk_py:$PATH"

  run gatk --java-options "$jopt" NVScoreVariants -R "$REF" -V "$hc" -O "$cnn" \
      --tensor-type reference
  (( DRY_RUN )) || gatk IndexFeatureFile -I "$cnn"

  # ---- 5. Tranche 过滤：只改 FILTER 列，不改位点 ----
  local hapmap
  hapmap=$(bundle_file hapmap_3.3.hg38 || true)
  [[ -n "$hapmap" ]] && ensure_vcf_index "$hapmap"
  if [[ -n "$hapmap" && -n "$known2" ]]; then
    run gatk --java-options "$jopt" FilterVariantTranches -V "$cnn" -O "$fin" \
        --info-key CNN_1D --snp-tranche 99.95 --indel-tranche 99.4 \
        --resource "$hapmap" --resource "$known2"
  else
    warn "缺少 hapmap_3.3.hg38，跳过 FilterVariantTranches（FILTER 未分层）"
    (( DRY_RUN )) || cp -f "$cnn" "$fin"
    (( DRY_RUN )) || cp -f "${cnn}.tbi" "${fin}.tbi" 2>/dev/null || true
  fi
  (( DRY_RUN )) || gatk IndexFeatureFile -I "$fin"

  log "GATK 完成: $fin"
}
