#!/usr/bin/env bash
# =============================================================================
#  src/32_stage_align.sh —— ③ bowtie2 比对
#
#  --local --very-fast-local：ATAC 的 reads 常在接头/低复杂度处截断，
#  允许局部比对能把这类 reads 救回来，整体比对率通常 >95%。
#  --no-mixed --no-discordant：只留正常配对，ATAC 不需要单端/跨染色体配对。
#  -X 2000：双端间距上限。设小了会把跨核小体的长片段丢掉。
#
#  输出 SAM 而不是直接管道进 BAM：③ 的产物就是 SAM（作业链上的"产出物"），
#  ④ 才做过滤和排序。磁盘紧张时把 KEEP_SAM 设 0，管道进 samtools。
# =============================================================================

stage_align() {
  section "③ bowtie2 比对"
  need bowtie2
  need samtools

  need_var "$BT2" "没找到 bowtie2 索引 —— 跑 deploy/deploy.sh bt2"
  # BT2 为空时上面已经提示过了，这里再拼 "$BT2.1.bt2" 只会打出个 ".1.bt2"
  if [[ -n "$BT2" ]]; then
    need_bt2 "$BT2"
  fi

  local sam="$MAP/${SAMPLE}.sam"
  local r1="$TRIM/${SAMPLE}.clean_R1.fastq.gz"
  local r2="$TRIM/${SAMPLE}.clean_R2.fastq.gz"
  need_file "$r1" "② 的清洗后 R1（先跑 --stages trim）"

  if (( KEEP_SAM )); then
    run bowtie2 "${BT2_OPTS[@]}" -p "$THREADS" -x "$BT2" -1 "$r1" -2 "$r2" \
        --rg-id "$SAMPLE" --rg "SM:$SAMPLE" --rg PL:ILLUMINA \
        -S "$sam" 2> "$MAP/${SAMPLE}.bowtie2.txt"
    (( DRY_RUN )) && { mark_done align; return 0; }

    local rate
    rate=$(awk '/overall alignment rate/ {gsub(/%/,"",$1); print $1; exit}' "$MAP/${SAMPLE}.bowtie2.txt")
    [[ -n "$rate" ]] && metric "03_bowtie2" "overall_alignment_rate" "$rate" "%" ">95"

    # chrM 占比：一遍 awk 同时拿到总比对数和 chrM 数
    # （不用 grep -c chrM：chrM_alt 也会被算进去）
    local counts total chrm
    counts=$(awk -F'\t' '!/^@/{t++; if($3=="chrM") m++} END{printf "%d\t%d", t, m+0}' "$sam")
    total=${counts%%$'\t'*}; chrm=${counts##*$'\t'}
    if [[ "${total:-0}" -gt 0 ]]; then
      local frac
      frac=$(awk -v m="$chrm" -v t="$total" 'BEGIN{printf "%.4f", m/t}')
      metric "03_bowtie2" "chrm_fraction" "$frac" "ratio" "<0.05"
      ok "比对率 ${rate}%，chrM 占 ${frac}（$chrm / $total）"
    fi
  else
    # 省磁盘路线：③④ 合并成一条管道，不落 SAM
    runpipe "bowtie2 ${BT2_OPTS[*]} -p $THREADS -x '$BT2' -1 '$r1' -2 '$r2' \
             --rg-id $SAMPLE --rg SM:$SAMPLE --rg PL:ILLUMINA \
             2> '$MAP/${SAMPLE}.bowtie2.txt' \
             | samtools view -@ $THREADS -b -q $MIN_MAPQ -e '$(contig_expr)' - \
             | samtools sort -@ $THREADS -m $SORT_MEM -T '$TMPD/.tmp_$SAMPLE' -o '$BAM_SRT' -"
    (( DRY_RUN )) && { mark_done align; return 0; }
    run samtools index -@ "$THREADS" "$BAM_SRT"
    local rate
    rate=$(awk '/overall alignment rate/ {gsub(/%/,"",$1); print $1; exit}' "$MAP/${SAMPLE}.bowtie2.txt")
    [[ -n "$rate" ]] && metric "03_bowtie2" "overall_alignment_rate" "$rate" "%" ">95"
    ok "比对率 ${rate}%（未落 SAM，③④ 合并执行）"
  fi

  mark_done align
}
