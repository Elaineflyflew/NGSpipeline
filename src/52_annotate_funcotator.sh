#!/usr/bin/env bash
# =============================================================================
#  src/52_annotate_funcotator.sh —— GATK Funcotator（对应 lab7）
#
#  同时出 VCF 和 MAF：VCF 给 IGV 和下游工具，MAF 直接丢进 Excel/表格看。
#  注意 Funcotator 要求参考的 contig 命名与数据源一致（都用 chr 前缀）。
# =============================================================================

annotate_funcotator() {
  local vcf="$1" tag="$2"
  need gatk
  require_arg "$FUNCOTATOR_DATA" "--funcotator-data" "Funcotator 数据源目录"
  checkdir "$FUNCOTATOR_DATA" "Funcotator 数据源"

  local jopt="-Xmx${JAVA_MEM}"
  local out="$ANN/${SAMPLE}.${tag}.funcotated.vcf.gz"
  local maf="$ANN/${SAMPLE}.${tag}.funcotated.maf"

  run gatk --java-options "$jopt" Funcotator \
      --variant "$vcf" --reference "$REF" --ref-version hg38 \
      --data-sources-path "$FUNCOTATOR_DATA" \
      --output "$out" --output-file-format VCF

  run gatk --java-options "$jopt" Funcotator \
      --variant "$vcf" --reference "$REF" --ref-version hg38 \
      --data-sources-path "$FUNCOTATOR_DATA" \
      --output "$maf" --output-file-format MAF

  (( DRY_RUN )) || tabix -f -p vcf "$out"
  log "Funcotator 完成: $out"
}
