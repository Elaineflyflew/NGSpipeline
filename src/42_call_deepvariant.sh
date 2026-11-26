#!/usr/bin/env bash
# =============================================================================
#  src/42_call_deepvariant.sh —— DeepVariant 路线（对应 lab8）
#
#  有 singularity 镜像就容器里跑，没有就用本机 run_deepvariant。
#  --model_type=WES 不能省：WGS 模型用在捕获数据上会有系统性偏差。
# =============================================================================

call_deepvariant() {
  section "⑤-C DeepVariant 变异检测"

  local d="$VAR/deepvariant"
  local vcf="$d/${SAMPLE}.dv.vcf.gz"
  local gvcf="$d/${SAMPLE}.dv.gvcf.gz"
  local region=()
  [[ -n "$BED" ]] && region=("--regions=$BED")

  if [[ -n "$DV_SIF" ]]; then
    (( DRY_RUN )) || [[ -s "$DV_SIF" ]] || die "singularity 镜像不存在: ${DV_SIF}（deploy/deploy.sh deepvariant）"
    local eng="singularity"
    command -v singularity >/dev/null 2>&1 || eng="apptainer"
    need "$eng"
    run "$eng" exec -B "$DV_BIND" "$DV_SIF" run_deepvariant \
        --reads="$MDUP" --ref="$REF" --model_type="$DV_MODEL" \
        --output_vcf="$vcf" --output_gvcf="$gvcf" \
        --num_shards="$THREADS" --vcf_stats_report=true \
        ${region[@]+"${region[@]}"}
  else
    need run_deepvariant
    run run_deepvariant --reads="$MDUP" --ref="$REF" --model_type="$DV_MODEL" \
        --output_vcf="$vcf" --output_gvcf="$gvcf" \
        --num_shards="$THREADS" --vcf_stats_report=true \
        ${region[@]+"${region[@]}"}
  fi

  log "DeepVariant 完成: $vcf"
}
