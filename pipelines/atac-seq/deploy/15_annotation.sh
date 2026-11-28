#!/usr/bin/env bash
# =============================================================================
#  deploy/15_annotation.sh —— 基因注释 GTF（⑪ 的 TSS 富集要用）
#
#  用 GENCODE 而不是 RefSeq：GENCODE 的 gene_type 字段干净，筛
#  protein_coding 只要一次字符串匹配；RefSeq 要处理一堆转录本别名。
# =============================================================================

: "${GENCODE_RELEASE:=43}"
: "${GENCODE_BASE:=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human}"

deploy_annotation() {
  section "基因注释 GTF（GENCODE v${GENCODE_RELEASE}）"
  check_space 2 "$ANN_DIR"

  local name="gencode.v${GENCODE_RELEASE}.annotation.gtf"
  local gz="$ANN_DIR/${name}.gz"
  local out="$ANN_DIR/${name}"

  dl "$GENCODE_BASE/release_${GENCODE_RELEASE}/${name}.gz" "$gz"

  if [[ -s "$out" ]] && (( ! FORCE )); then
    log "已存在，跳过解压: $out"
  elif (( DRY_RUN )); then
    printf '[dry-run] gunzip -c %s > %s\n' "$gz" "$out"
  else
    log "解压"
    gunzip -c "$gz" > "$out" || die "解压失败: $gz"
    rm -f "$gz"
  fi

  ok "GTF: $out（运行时按 *.gtf 自动探测）"
}
