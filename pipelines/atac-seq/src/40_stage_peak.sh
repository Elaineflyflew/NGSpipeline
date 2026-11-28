#!/usr/bin/env bash
# =============================================================================
#  src/40_stage_peak.sh —— ⑧⑨ peak calling + 去黑名单
#
#  三条路线各 call 一次：
#    · 全部片段（srt2）  常规 peak，报告主结果用这份
#    · NFR              只用裸 DNA 片段，信噪比最高，峰更"尖"
#    · mono             跨单个核小体的片段，反映核小体定位
#  三份 peak 的数目差异本身就是结论：NFR 应该比 mono 多，
#  如果反过来，说明片段分档或偏移校正有问题。
#
#  MACS2 参数说明：
#   -f BAMPE       用双端 reads 的真实跨度定义片段，不要再用 --extsize 猜
#   --nomodel      BAMPE 模式下 MACS 本来就会自动关掉建模，显式写出来是为了读代码的人
#   -g hs          人类有效基因组大小（2.7e9）
#   不加 -c        ATAC 没有 IP/input 对照，加 -c 反而引入噪声
#   不传 --shift   BAMPE 模式下 MACS 会把 shift 强制归零（源码 OptValidator.py），
#                  传了也不生效还不报错 —— 偏移校正要在 ⑤ 的 BAM 层面做完
# =============================================================================

# peak 覆盖的基因组比例（snakePipes 的 MACS2_QC 里有一项就是这个）
peak_genome_pct() {
  local peaks="$1" gsize="$2"
  awk -v g="$gsize" '{ s += $3 - $2 } END { printf "%.2f", (g > 0 ? 100 * s / g : 0) }' "$peaks"
}

genome_size() {
  local fai="${FASTA}.fai"
  if [[ -s "$fai" ]]; then
    awk '{ s += $2 } END { printf "%d", s }' "$fai"
  else
    printf '3100000000'   # hg38 主染色体的近似总长，仅用于估算百分比
  fi
}

call_peaks_one() {
  local bam="$1" name="$2" out="$3"

  run macs2 callpeak \
      -t "$bam" \
      -f "$PEAK_FORMAT" --nomodel \
      -g "$EFFECTIVE_GENOME" \
      -q "$MACS2_QVALUE" \
      -n "$name" \
      --outdir "$PEAK"

  # 去 ENCODE 黑名单：-v 是反选，保留与黑名单不重叠的 peak
  runpipe "bedtools intersect -v -a '$PEAK/${name}_peaks.narrowPeak' -b '$BLACKLIST' > '$out'"
}

stage_peak() {
  section "⑧ MACS2 callpeak（全部 / NFR / mono）"
  need macs2; need bedtools

  checkfile "$SRT2" "⑥ 产出的 srt2.bam"
  checkfile "$NFR"  "⑦ 产出的 NFR BAM"
  checkfile "$MONO" "⑦ 产出的 mono BAM"
  checkfile "$BLACKLIST" "ENCODE blacklist"

  call_peaks_one "$SRT2" "$SAMPLE"      "$PEAKS"
  call_peaks_one "$NFR"  "${SAMPLE}_NFR"  "$PEAKS_NFR"
  call_peaks_one "$MONO" "${SAMPLE}_mono" "$PEAKS_MONO"

  section "⑨ 去黑名单 + peak 统计"
  if (( ! DRY_RUN )); then
    local gsize
    gsize=$(genome_size)
    {
      printf 'set\tpeaks\tfiltered\tgenome_pct\n'
      local tag raw filt
      for tag in all nfr mono; do
        case "$tag" in
          all)  raw="$PEAKS_RAW";   filt="$PEAKS" ;;
          nfr)  raw="$PEAK/${SAMPLE}_NFR_peaks.narrowPeak";   filt="$PEAKS_NFR" ;;
          mono) raw="$PEAK/${SAMPLE}_mono_peaks.narrowPeak";  filt="$PEAKS_MONO" ;;
        esac
        printf '%s\t%s\t%s\t%s\n' "$tag" \
          "$(wc -l < "$raw"  | tr -d ' ')" \
          "$(wc -l < "$filt" | tr -d ' ')" \
          "$(peak_genome_pct "$filt" "$gsize")"
      done
    } > "$PEAK/peak_counts.tsv"
    log "peak 数: $(awk -F'\t' 'NR>1{printf "%s=%s ", $1, $3}' "$PEAK/peak_counts.tsv")"
  fi

  mark_done peak
}
