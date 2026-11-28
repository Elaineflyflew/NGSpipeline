#!/usr/bin/env bash
# =============================================================================
#  deploy/30_homer.sh —— HOMER + hg38 数据
#
#  HOMER 不是"装个二进制"就完事，它是 perl 脚本 + 一套数据目录：
#    1. 下 configureHomer.pl
#    2. perl configureHomer.pl -install            # 装 HOMER 本体到 ./homer
#    3. perl configureHomer.pl -install hg38       # 下基因组序列（~3 G）
#  第 3 步不做的话，annotatePeaks.pl 会报 "Could not find genome hg38"。
#
#  known.motifs 是随本体一起下的，路径：homer/data/knownTFs/known.motifs
# =============================================================================

: "${HOMER_BASE:=http://homer.ucsd.edu/homer}"

deploy_homer() {
  section "HOMER + hg38"
  need perl

  mkdir -p "$HOMER_DIR"
  local cfg="$HOMER_DIR/configureHomer.pl"

  dl "$HOMER_BASE/configureHomer.pl" "$cfg"

  if (( DRY_RUN )); then
    printf '[dry-run] (cd %s && perl configureHomer.pl -install)\n' "$HOMER_DIR"
    printf '[dry-run] (cd %s && perl configureHomer.pl -install hg38)\n' "$HOMER_DIR"
    printf '[dry-run] known.motifs 会落在 %s/data/knownTFs/known.motifs\n' "$HOMER_DIR"
    return 0
  fi

  # -install 不带参数：只装本体，不下基因组（避免一次拉几十 G）
  if [[ -x "$HOMER_DIR/bin/annotatePeaks.pl" ]]; then
    log "HOMER 本体已存在，跳过"
  else
    log "安装 HOMER 本体"
    ( cd "$HOMER_DIR" && perl configureHomer.pl -install ) || die "HOMER 本体安装失败"
  fi

  if [[ -d "$HOMER_DIR/data/genomes/hg38" ]]; then
    log "hg38 数据已存在，跳过"
  else
    log "下载 hg38 基因组数据（约 3 G）"
    ( cd "$HOMER_DIR" && perl configureHomer.pl -install hg38 ) || die "hg38 数据下载失败"
  fi

  # 把 HOMER 的 bin 挂到 extratools/bin，这样 src/ 的 find_tool 不用改配置就能找到
  mkdir -p "$BIN_DIR"
  local t
  for t in annotatePeaks.pl findMotifsGenome.pl; do
    [[ -x "$HOMER_DIR/bin/$t" ]] && ln -sf "$HOMER_DIR/bin/$t" "$BIN_DIR/$t"
  done

  ok "HOMER: $HOMER_DIR"
  ok "known.motifs: $HOMER_DIR/data/knownTFs/known.motifs"
}
