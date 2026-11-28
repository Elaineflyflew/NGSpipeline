#!/usr/bin/env bash
# =============================================================================
#  deploy/50_igv.sh —— IGV（⑮ 截图）
#
#  ⑮ 有两条路：
#    有 igv        → 跑 batch 脚本，xvfb-run 包一层就能在无显示器的服务器上出图
#    没有 igv      → 退回 pyGenomeTracks（conda 装，更轻）
#  两条都没有也行：batch 文件照样生成，拷到本地双击就能截图。
# =============================================================================

: "${IGV_VERSION:=2.19.5}"
: "${IGV_BASE:=https://data.broadinstitute.org/igv/projects/downloads}"

deploy_igv() {
  section "IGV ${IGV_VERSION}（可选）"

  local os url zip
  case "$(uname -s)" in
    Linux)  os="Linux"; url="$IGV_BASE/2.19/IGV_Linux_${IGV_VERSION}.zip" ;;
    Darwin) os="Mac";   url="$IGV_BASE/2.19/IGV_${IGV_VERSION}.zip" ;;
    *)      warn "不认识的系统 $(uname -s)，只装 pyGenomeTracks"; os="" ;;
  esac

  if [[ -n "$os" ]]; then
    zip="$IGV_DIR/IGV_${os}_${IGV_VERSION}.zip"
    mkdir -p "$IGV_DIR"
    dl "$url" "$zip"
    if (( DRY_RUN )); then
      printf '[dry-run] unzip -q %s -d %s\n' "$zip" "$IGV_DIR"
    else
      unzip -q -o "$zip" -d "$IGV_DIR" || die "解压失败: $zip"
      # 解压出来是 IGV_<ver>/igv.sh，软链到 extratools/bin/igv
      local sh
      sh=$(find "$IGV_DIR" -maxdepth 2 -name 'igv.sh' -o -maxdepth 2 -name 'igv' -type f 2>/dev/null | head -1)
      if [[ -n "$sh" ]]; then
        chmod +x "$sh"
        mkdir -p "$BIN_DIR"
        ln -sf "$sh" "$BIN_DIR/igv"
        ok "igv → $BIN_DIR/igv"
      fi
    fi
  fi

  # 服务器上没有 X11，xvfb-run 是让 IGV 能出图的唯一办法
  if ! command -v xvfb-run >/dev/null 2>&1; then
    warn "没有 xvfb-run：服务器上跑 IGV 需要它（apt install xvfb / conda install -c conda-forge xorg-xvfb）"
  fi
  if ! find_tool pyGenomeTracks >/dev/null 2>&1; then
    log "没有 pyGenomeTracks（可选兜底）：conda install -c bioconda pygenometracks"
  fi
  ok "IGV 目录: $IGV_DIR"
}
