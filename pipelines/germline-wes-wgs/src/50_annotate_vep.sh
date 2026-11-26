#!/usr/bin/env bash
# =============================================================================
#  src/50_annotate_vep.sh —— VEP 注释（对应 lab6）
#
#  cache 路径有两种常见写法，这里统一规整：
#    .../vep                        → 正确（其下有 homo_sapiens/）
#    .../vep/homo_sapiens/115_GRCh38 → 多了一层，自动往上剥
# =============================================================================

annotate_vep() {
  local vcf="$1" tag="$2"
  need vep; need bgzip; need tabix

  local cache="$VEP_CACHE"
  [[ "$cache" == *"/homo_sapiens/"* ]] && cache="${cache%%/homo_sapiens/*}"
  checkdir "$cache/homo_sapiens" "VEP cache（其下应有 homo_sapiens/）"

  local out="$ANN/${SAMPLE}.${tag}.vep.vcf"
  run vep --fork "$THREADS" --cache --offline --dir_cache "$cache" \
      --assembly "$VEP_ASSEMBLY" --species homo_sapiens \
      -i "$vcf" -o "$out" --vcf --force_overwrite \
      --symbol --canonical --biotype --hgvs --numbers --variant_class \
      --sift p --polyphen p --af_gnomade --max_af

  if (( ! DRY_RUN )); then
    bgzip -f "$out"
    tabix -f -p vcf "${out}.gz"
  fi
  log "VEP 完成: ${out}.gz"
}
