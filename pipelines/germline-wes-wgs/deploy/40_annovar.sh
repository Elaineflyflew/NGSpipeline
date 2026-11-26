#!/usr/bin/env bash
# =============================================================================
#  deploy/40_annovar.sh —— ANNOVAR + humandb
#
#  ANNOVAR 需要学术注册才能下载，脚本不能替你点同意，所以分两步：
#    1. 你把 annovar 压缩包放进 extratools/ 或 --tarball 指定路径
#    2. 本脚本解压 + 用 annotate_variation.pl 下载 humandb
#
#  没给压缩包时只打印指引，不报错（避免 deploy all 被卡住）。
# =============================================================================
set -euo pipefail

# 与 config/pipeline.conf 的 ANNOVAR_PROTOCOLS 保持一致
ANNOVAR_DBS="${ANNOVAR_DBS:-refGene cytoBand exac03 gnomad30_genome avsnp150 dbnsfp42a}"

deploy_annovar() {
  local tarball=""
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --tarball) tarball="$2"; shift 2 ;;
      -h|--help) echo "用法: $0 [--tarball annovar.tar.gz]"; return 0 ;;
      *) die "未知参数: $1" ;;
    esac
  done

  mkdir -p "$ANNOVAR_DIR"

  # ---- 1. 解压 ----
  if [[ ! -f "$ANNOVAR_DIR/annotate_variation.pl" ]]; then
    if [[ -z "$tarball" ]]; then
      # 本地已经有下载好的包就自动用
      tarball=$(first_existing "$TOOLS_DIR"/annovar*.tar.gz "$TOOLS_DIR"/annovar*.tar.bz2 || true)
    fi
    if [[ -z "$tarball" || ! -s "$tarball" ]]; then
      cat <<'EOF'

  ANNOVAR 需要注册后手动下载（不能脚本代下）：
    1. 打开 https://www.openbioinformatics.org/annovar/annovar_download_form.php
    2. 填表后邮件会收到下载链接，下载 annovar.latest.tar.gz
    3. 放到 extratools/ 下，然后重跑：
         bash deploy/40_annovar.sh --tarball extratools/annovar.latest.tar.gz

  注意：没有 ANNOVAR 不影响其他步骤，--annotator 用 vep 或 funcotator 即可。

EOF
      return 0
    fi
    log "解压 $tarball"
    tar -xf "$tarball" -C "$ANNOVAR_DIR" --strip-components=1
  else
    log "ANNOVAR 已就绪: $ANNOVAR_DIR"
  fi

  # ---- 2. humandb ----
  local perl="$ANNOVAR_DIR/annotate_variation.pl"
  [[ -f "$perl" ]] || die "解压后找不到 annotate_variation.pl，检查压缩包结构"
  mkdir -p "$ANNOVAR_DIR/humandb"

  local db
  for db in $ANNOVAR_DBS; do
    if (( DRY_RUN )); then
      printf '[dry-run] perl %s -buildver hg38 -downdb -webfrom annovar %s humandb/\n' "$perl" "$db"
      continue
    fi
    # dbnsfp 系列单库就 10G+，先提醒
    case "$db" in
      dbnsfp*) check_space 15 "$ANNOVAR_DIR/humandb" ;;
    esac
    log "下载 humandb: $db"
    perl "$perl" -buildver hg38 -downdb -webfrom annovar "$db" "$ANNOVAR_DIR/humandb" \
      || warn "$db 下载失败（网络或该库已改名），可稍后单独重跑本脚本"
  done

  log "ANNOVAR 就绪: $ANNOVAR_DIR/humandb"
}

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
  source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
  deploy_annovar "$@"
fi
