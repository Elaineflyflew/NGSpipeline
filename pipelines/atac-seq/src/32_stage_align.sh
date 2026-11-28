#!/usr/bin/env bash
# =============================================================================
#  src/32_stage_align.sh —— ③ bowtie2 比对 → SAM
#
#  ATAC-seq 为什么用 --local 而不是 --end-to-end：
#  开放区的 read 常跨越 Tn5 切口或小 indel，末端允许错配能显著提高比对率。
#
#  --no-mixed / --no-discordant：ATAC 只关心正常配对的片段，
#  单端比对和跨染色体配对大多是噪声，留着只会干扰后续 TLEN 分档。
#
#  注意 --local 已经隐含在 --very-fast-local 里，两个一起写是冗余的；
#  保留是为了和课程脚本的命令行逐字对得上，方便对照。
# =============================================================================

stage_align() {
  section "③ 比对 bowtie2"
  need bowtie2

  # 有修剪后的就用，没有就回退到原始 FASTQ（--stages 里跳过 trim 时走这里）
  local in1="$R1" in2="$R2"
  local c1="$CLEAN/${SAMPLE}_R1.clean.fq.gz"
  local c2="$CLEAN/${SAMPLE}_R2.clean.fq.gz"
  if [[ -s "$c1" && -s "$c2" ]]; then in1="$c1"; in2="$c2"; fi

  checkfile "$BT2_INDEX.1.bt2" "bowtie2 索引"

  runpipe "bowtie2 --local $BOWTIE2_PRESET \
      --no-mixed --no-discordant --phred33 \
      -p $THREADS -x '$BT2_INDEX' \
      -1 '$in1' -2 '$in2' \
      -S '$SAM' \
      2> '$MAP/${SAMPLE}_bowtie2.log'"

  if (( ! DRY_RUN )); then
    [[ -s "$SAM" ]] || die "bowtie2 没有产出 SAM: $SAM"
    local rate
    rate=$(grep -oE '[0-9.]+% overall alignment rate' "$MAP/${SAMPLE}_bowtie2.log" \
           | grep -oE '^[0-9.]+' | head -1)
    printf 'overall_alignment_rate\t%s\n' "${rate:--}" \
      > "$MAP/${SAMPLE}.align_rate.txt"
    log "比对率: ${rate:--}%"
  fi

  mark_done align
}
