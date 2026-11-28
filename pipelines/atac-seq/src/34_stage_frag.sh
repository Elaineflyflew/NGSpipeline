#!/usr/bin/env bash
# =============================================================================
#  src/34_stage_frag.sh —— ⑦ 按 |TLEN| 分箱
#
#  三档：NFR <147 bp / 单核小体 147–294 / 多核小体 >294
#  147 bp 是"缠绕一个核小体"的长度，所以这个阈值只在**已做 Tn5 偏移**的
#  srt2.bam 上才成立。用未偏移的 srt.bam 分档，判定窗口整体错位约 9 bp，
#  NFR 档会系统性偏少——这是最容易被忽略的一个坑。
#
#  产物：三份 BAM（喂 ⑧ 的三条 peak calling 路线）+ 片段分布图（给人看）
# =============================================================================

# 各档的参考区间。多核小体档没有硬阈值：它的用途是看周期性，不是判合格。
frag_bench() {
  case "$1" in
    nfr)  printf '>0.30' ;;
    mono) printf '0.20-0.40' ;;
    *)    printf '' ;;
  esac
}

stage_frag() {
  section "⑦ TLEN 分箱（NFR / 单核小体 / 多核小体）"
  need samtools
  need_file "$BAM_SRT2" "⑥ 的 srt2.bam（先跑 --stages bam）"

  run "$PYTHON" "$PROJECT_ROOT/src/py/frag_bin.py" \
      --bam "$BAM_SRT2" --outdir "$FRAG" --bamdir "$BAM" \
      --sample "$SAMPLE" \
      --nfr-max "$TLEN_NFR_MAX" --mono-max "$TLEN_MONO_MAX" \
      --threads "$THREADS"
  (( DRY_RUN )) && { mark_done frag; return 0; }

  # fragbin_stats.tsv: bin / fragments / fraction / mean_length / median_length
  local stats="$FRAG/${SAMPLE}.fragbin_stats.tsv"
  local bin frags frac mean med
  while IFS=$'\t' read -r bin frags frac mean med; do
    [[ "$bin" == "bin" || -z "$bin" ]] && continue
    metric "07_fragbin" "${bin}_fragments" "$frags" "fragments" ""
    metric "07_fragbin" "${bin}_fraction"  "$frac"  "ratio"     "$(frag_bench "$bin")"
    metric "07_fragbin" "${bin}_median_length" "$med" "bp" ""
  done < "$stats"

  local nfr mono poly
  nfr=$(awk -F'\t' '$1=="nfr"{print $3}'  "$stats")
  mono=$(awk -F'\t' '$1=="mono"{print $3}' "$stats")
  poly=$(awk -F'\t' '$1=="poly"{print $3}' "$stats")
  ok "NFR ${nfr} / 单核小体 ${mono} / 多核小体 ${poly}"
  mark_done frag
}
