#!/usr/bin/env bash
# =============================================================================
#  src/40_call_bcftools.sh —— bcftools 路线（对应 lab6）
#
#  命令链：mpileup → call -mv → filter 'QUAL>10 && DP>10' → 拆 SNP/indel
#  -d 1000 是关键：默认 250 会削掉 WES 目标区的高深度位点，AF 会失真。
# =============================================================================

call_bcftools() {
  section "⑤-A bcftools 变异检测"
  need bcftools

  local d="$VAR/bcftools"
  local mpi="$d/${SAMPLE}.bcftools.mpileup.vcf.gz"
  local raw="$d/${SAMPLE}.bcftools.raw.vcf.gz"
  local fin="$d/${SAMPLE}.bcftools.final.vcf.gz"
  local bedopt=()
  [[ -n "$BED" ]] && bedopt=(-R "$BED")

  run bcftools mpileup --threads "$THREADS" -f "$REF" -q "$MIN_MAPQ" -Q 20 -d 1000 \
      -a FORMAT/AD,FORMAT/DP,INFO/AD ${bedopt[@]+"${bedopt[@]}"} \
      -Oz -o "$mpi" "$MDUP"
  run bcftools index -t -f "$mpi"

  run bcftools call -mv --threads "$THREADS" -Oz -o "$raw" "$mpi"
  run bcftools index -t -f "$raw"

  run bcftools filter --threads "$THREADS" \
      -i "QUAL>${MIN_QUAL} && DP>${MIN_DP}" -Oz -o "$fin" "$raw"
  run bcftools index -t -f "$fin"

  run bcftools view -v snps   -Oz -o "$d/${SAMPLE}.bcftools.snps.vcf.gz"   "$fin"
  run bcftools view -v indels -Oz -o "$d/${SAMPLE}.bcftools.indels.vcf.gz" "$fin"

  if (( ! DRY_RUN )); then
    bcftools index -t -f "$d/${SAMPLE}.bcftools.snps.vcf.gz"
    bcftools index -t -f "$d/${SAMPLE}.bcftools.indels.vcf.gz"
    bcftools stats "$fin" > "$d/${SAMPLE}.bcftools.stats.txt"
    if command -v plot-vcfstats >/dev/null 2>&1; then
      plot-vcfstats -p "$d/plot" "$d/${SAMPLE}.bcftools.stats.txt" >/dev/null 2>&1 || true
    fi
  fi
  log "bcftools 完成: $fin"
}
