#!/usr/bin/env bash
# =============================================================================
#  src/31_stage_trim.sh —— ② 接头修剪 + 质量过滤，然后复检
#
#  ATAC 的两个注意点：
#    · 不开 --trim_poly_g：非 NovaSeq 建库没有 polyG 偏好，开了会误削真实碱基
#    · -q 不压到 20 以上：Tn5 插入位点附近的低质量碱基被削掉，切割位点会偏移
#  修剪完必须再看一遍 FastQC，否则"接头掉了"只是假设，不是证据。
# =============================================================================

stage_trim() {
  section "② fastp 修剪 + 清洗后复检"
  need fastp
  need fastqc
  fastqc_probe

  local out1="$TRIM/${SAMPLE}.clean_R1.fastq.gz"
  local out2="$TRIM/${SAMPLE}.clean_R2.fastq.gz"
  local json="$TRIM/${SAMPLE}.fastp.json"
  local html="$TRIM/${SAMPLE}.fastp.html"

  run fastp -i "$R1" -I "$R2" -o "$out1" -O "$out2" \
            "${FASTP_OPTS[@]}" -w "$THREADS" -j "$json" -h "$html"

  run fastqc -t "$THREADS" ${FASTQC_EXTRA[@]+"${FASTQC_EXTRA[@]}"} --outdir "$QCC" "$out1" "$out2"

  # fastp JSON → 五列指标（02_fastp.tsv）
  runpipe "'$PYTHON' '$PROJECT_ROOT/src/py/fastp_stats.py' --json '$json' --out '$METRIC/02_fastp.tsv' --step 02_fastp"

  (( DRY_RUN )) && { mark_done trim; return 0; }

  local kept q30
  kept=$(awk -F'\t' '$2=="reads_retained"{print $3}' "$METRIC/02_fastp.tsv" 2>/dev/null | head -1)
  q30=$(awk  -F'\t' '$2=="q30_after_trim"{print $3}'  "$METRIC/02_fastp.tsv" 2>/dev/null | head -1)
  ok "保留率 ${kept:-?}，Q30（清洗后）${q30:-?}"
  mark_done trim
}
