#!/usr/bin/env bash
# =============================================================================
#  src/33_stage_bam.sh —— ④ 过滤排序  ⑤ Tn5 偏移  ⑥ 重排
#
#  三步为什么必须按这个顺序：
#    ④ 先去掉 chrM 和低 MAPQ。chrM 占 20–30% 是 ATAC 常态，不去掉后面
#      的 TLEN 分箱和 peak calling 都会被线粒体信号淹没。
#    ⑤ alignmentSieve --ATACshift 把每条 read 的 5' 端 +4 / 3' 端 −5，
#      还原 Tn5 的真实切割位点。不做这一步，peak 会整体偏约 9 bp。
#    ⑥ shift 之后**必须重排**：alignmentSieve 只保留前 11 个字段（tag 全丢），
#      而且坐标被挪动后顺序已乱，直接 index 会得到错误的区间索引。
# =============================================================================

stage_bam() {
  section "④⑤⑥ BAM 过滤 → Tn5 偏移 → 重排"
  need samtools
  need alignmentSieve

  local expr; expr=$(contig_expr)

  # ---- ④ 过滤 + 排序 ----
  if (( KEEP_SAM )); then
    local sam="$MAP/${SAMPLE}.sam"
    need_file "$sam" "③ 的 SAM（先跑 --stages align）"
    runpipe "samtools view -@ $THREADS -b -q $MIN_MAPQ -e '$expr' '$sam' \
             | samtools sort -@ $THREADS -m $SORT_MEM -T '$TMPD/.tmp_$SAMPLE' -o '$BAM_SRT' -"
  else
    log "③ 已把过滤合并进比对管道，跳过重复过滤"
  fi
  run samtools index -@ "$THREADS" "$BAM_SRT"

  # ---- ④ 指标（要读文件，dry-run 跳过）----
  local total="" dup_pct=""
  if (( ! DRY_RUN )); then
    samtools flagstat -@ "$THREADS" "$BAM_SRT" > "$BAM/${SAMPLE}.srt.flagstat.txt"
    samtools stats    -@ "$THREADS" "$BAM_SRT" > "$BAM/${SAMPLE}.srt.stats.txt"

    local fstat="$BAM/${SAMPLE}.srt.flagstat.txt"
    local pp dup_n ins
    # flagstat 的行形如
    #   "20000000 + 0 in total (QC-passed reads + QC-failed reads)"
    #   "25177548 + 0 properly paired (25177548 : 100.00%)"
    # 注意 properly paired 那行的百分号在最后一个字段，而且老版本 samtools 写成
    # "(100.00% : N/A)" —— 所以不能取固定列号（取 $5 会拿到单词 "paired"，
    # 之前就写出过 proper_pair_rate=paired 这种值）。改成扫第一个含 % 的字段。
    total=$(awk '$2=="+" && $5=="total" {print $1; exit}' "$fstat")
    pp=$(awk '/properly paired/ {
                for (i = 1; i <= NF; i++)
                  if ($i ~ /%/) { gsub(/[()%]/, "", $i); print $i; exit }
              }' "$fstat")
    dup_n=$(awk '$4=="duplicates" {print $1; exit}' "$fstat")
    dup_pct=$(awk -v d="${dup_n:-0}" -v t="${total:-0}" 'BEGIN{printf "%.2f", (t>0? d*100/t : 0)}')
    ins=$(stats_field "$BAM/${SAMPLE}.srt.stats.txt" "insert size average")

    metric "04_bam_filter" "reads_after_filter" "${total:-NA}" "reads" ""
    metric "04_bam_filter" "proper_pair_rate"   "${pp:-NA}"    "%"     ">80"
    # 本流程**不做 markdup**（和课程脚本一致；⑧ MACS2 用 --keep-dup all）。
    # 没跑 markdup 时 bowtie2 不会给任何 read 打 0x400，flagstat 的
    # "duplicates" 恒为 0，报成 0% 会被读成"没有重复"——这是假信息。
    # 所以这里显式写 NA，真实重复水平看 ② 的 fastp 估计值（dup 22.8%）。
    metric "04_bam_filter" "duplicates_flagged" "$dup_pct"    "%"     ""
    metric "04_bam_filter" "duplicate_rate"     "NA"          "%"     ""
    metric "04_bam_filter" "insert_size_mean"   "${ins:-NA}"   "bp"    "150-500"
  fi

  # ---- ⑤ Tn5 偏移 ----
  # deeptools 每个线程开一个临时 BAM，线程数超过文件句柄上限就会
  # "Too many open files"（报的还是 /tmp 里的文件名，很容易误判成磁盘问题）。
  local st="$SIEVE_THREADS"
  (( st > THREADS )) && st="$THREADS"
  run alignmentSieve -p "$st" --ATACshift -b "$BAM_SRT" -o "$BAM_SHIFT"

  # ---- ⑥ 重排 ----
  run samtools sort -@ "$THREADS" -m "$SORT_MEM" -T "$TMPD/.tmp2_$SAMPLE" \
        -o "$BAM_SRT2" "$BAM_SHIFT"
  run samtools index -@ "$THREADS" "$BAM_SRT2"

  (( DRY_RUN )) && { mark_done bam; return 0; }

  # 字段数必须是 11：少一个说明 tag 被丢、后续 bedtools/HOMER 读不出链方向
  # 不能用 `samtools view | head -1`：pipefail 下 head 提前退出会让 samtools
  # 收到 SIGPIPE，整条命令判失败
  local nf=""
  samtools view "$BAM_SRT2" 2>/dev/null | head -1 > "$TMPD/.f1.sam" || true
  nf=$(awk -F'\t' 'NF{print NF; exit}' "$TMPD/.f1.sam")
  rm -f "$TMPD/.f1.sam"
  [[ "${nf:-0}" == "11" ]] || die "⑥ 之后 SAM 字段数为 ${nf}（应为 11）—— alignmentSieve 的输出被后续工具改过"
  metric "06_resort" "bam_field_count" "11" "fields" "11-11"

  ok "④⑤⑥ 完成：$(basename "$BAM_SRT2")（reads ${total}，重复率 ${dup_pct}%）"
  mark_done bam
}
