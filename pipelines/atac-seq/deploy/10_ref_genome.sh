#!/usr/bin/env bash
# =============================================================================
#  deploy/10_ref_genome.sh —— 参考基因组 fasta + bowtie2 索引
#
#  为什么分两步（genome / bt2）：fasta 是 1 G，建索引要再花 1 G 和 1 小时。
#  换比对器（bwa/hisat2）时只需要重跑索引那一步，不用重下 fasta。
# =============================================================================

: "${UCSC_BASE:=https://hgdownload.soe.ucsc.edu/goldenPath/hg38}"
: "${GENOME_NAME:=hg38}"

deploy_genome() {
  section "参考基因组 fasta（${GENOME_NAME}）"
  check_space 4 "$GENOME_DIR"

  local gz="$GENOME_DIR/${GENOME_NAME}.fa.gz"
  local fa="$GENOME_DIR/${GENOME_NAME}.fa"

  dl "$UCSC_BASE/bigZips/${GENOME_NAME}.fa.gz" "$gz"

  if [[ -s "$fa" ]] && (( ! FORCE )); then
    log "已存在，跳过解压: $fa"
  elif (( DRY_RUN )); then
    printf '[dry-run] gunzip -c %s > %s\n' "$gz" "$fa"
  else
    log "解压（约 3 G，几分钟）"
    gunzip -c "$gz" > "$fa" || die "解压失败: $gz"
    rm -f "$gz"
  fi

  need samtools
  run samtools faidx "$fa"
  ok "fasta: $fa"
}

deploy_bt2() {
  section "bowtie2 索引"
  need bowtie2-build

  local fa; fa=$(find_genome)
  [[ -n "$fa" ]] || die "先跑 deploy/deploy.sh genome"
  check_space 5 "$BT2_DIR"

  mkdir -p "$BT2_DIR"
  run bowtie2-build --threads "${THREADS:-8}" "$fa" "$BT2_DIR/$GENOME_NAME"
  ok "索引前缀: $BT2_DIR/$GENOME_NAME（运行时自动探测，不用改配置）"
}
