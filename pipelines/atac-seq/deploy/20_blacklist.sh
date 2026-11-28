#!/usr/bin/env bash
# =============================================================================
#  deploy/20_blacklist.sh —— ENCODE hg38 blacklist v2
#
#  只有 1 MB，但 ⑨ 和 ⑪ 都要用：
#    ⑨ 去掉落在黑名单里的 peak（着丝粒/端粒/rDNA 高拷贝区）
#    ⑪ bamCoverage --blackListFileName 让信号图不出现假高峰
# =============================================================================

: "${BLACKLIST_URL:=https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz}"

deploy_blacklist() {
  section "ENCODE blacklist v2"
  mkdir -p "$BL_DIR"

  local gz="$BL_DIR/hg38-blacklist.v2.bed.gz"
  local out="$BL_DIR/hg38-blacklist.v2.bed"

  if [[ -s "$out" ]] && (( ! FORCE )); then
    log "已存在，跳过: $out"
    return 0
  fi

  dl "$BLACKLIST_URL" "$gz"

  if (( DRY_RUN )); then
    printf '[dry-run] gunzip -c %s > %s\n' "$gz" "$out"
    return 0
  fi
  gunzip -c "$gz" > "$out" || die "解压失败: $gz"
  rm -f "$gz"
  ok "黑名单: $out（$(wc -l < "$out") 个区间）"
}
