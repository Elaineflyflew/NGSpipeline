#!/usr/bin/env bash
# =============================================================================
#  src/53_stage_annotate.sh —— ⑥ 注释调度
#
#  caller × annotator 的笛卡尔积，每个组合单独一个完成标记。
#  某个注释库没配好只会跳过那一格，不影响其他组合，也不会把前面的比对结果废掉。
# =============================================================================

stage_annotate() {
  if (( ${#ANNOTATORS[@]} == 0 )); then
    log "⑥ 跳过注释（--annotator none）"
    mark_done annotate
    return 0
  fi

  local c a vcf key
  for c in "${CALLERS[@]}"; do
    vcf=$(final_vcf "$c")
    if (( ! DRY_RUN )) && [[ ! -s "$vcf" ]]; then
      warn "$c 的 VCF 不存在，跳过其注释: $vcf"
      continue
    fi
    for a in "${ANNOTATORS[@]}"; do
      key="annotate_${a}_${c}"
      if stage_done "$key"; then
        log "⑥ [$a/$c] 已完成，跳过"
        continue
      fi
      case "$a" in
        vep)        annotate_vep        "$vcf" "$c" ;;
        annovar)    annotate_annovar    "$vcf" "$c" ;;
        funcotator) annotate_funcotator "$vcf" "$c" ;;
        *) die "内部错误：未知 annotator $a" ;;
      esac
      mark_done "$key"
    done
  done
  mark_done annotate
}
