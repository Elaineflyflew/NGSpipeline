#!/usr/bin/env bash
# =============================================================================
#  src/43_stage_call.sh —— ⑤ 变异检测调度
#
#  只负责"按顺序调谁"，具体命令在 40/41/42 里。
#  每个 caller 单独落一个 call_<name>.done，所以某一条路线失败后
#  重跑只会补那一条，不会把已经跑完的另外两条重算。
# =============================================================================

stage_call() {
  local c
  for c in "${CALLERS[@]}"; do
    if stage_done "call_$c"; then
      log "⑤ [$c] 已完成，跳过"
      continue
    fi
    case "$c" in
      bcftools)    call_bcftools ;;
      gatk)        call_gatk ;;
      deepvariant) call_deepvariant ;;
      *) die "内部错误：未知 caller $c" ;;
    esac
    mark_done "call_$c"
  done
  mark_done call
}
