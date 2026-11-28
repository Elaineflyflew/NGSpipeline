#!/usr/bin/env bash
# =============================================================================
#  src/51_stage_motif.sh —— ⑬ motif 富集（只比对已知库）
#
#  为什么先取 top N：全量 7 万 peak 跑 known motif 要 2 小时，取信号最强的
#  1 万既省时间，又不改变富集结论——motif 富集本来就是被最强的那批 peak
#  驱动的。
#
#  为什么加 -mknown：只跑"已知库比对"这一半，结果就是 knownResults.txt。
#  想连 de novo discovery 一起跑就设 MOTIF_DENOVO=1（不加 -nomotif）。
#
#  实测提醒：不同 HOMER 版本对 -nomotif 的响应不一致——本仓库在 HOMER 4.11 上
#  跑时，发现步骤仍然执行，只是省掉了最后那轮 de novo 富集统计。好消息是
#  1 万个 peak 的 de novo 只要 3 分钟左右（7 万 peak 才要 1–2 小时），
#  所以这个开关留着不亏，但也别指望它一定省下全部时间。
# =============================================================================

stage_motif() {
  section "⑬ HOMER known motif 富集"
  need findMotifsGenome.pl
  need_var "$HOMER_MOTIFS" "没找到 HOMER 的 known.motifs —— 跑 deploy/deploy.sh homer，或用 --homer-motifs 指定"

  # de novo 开关：默认关（-nomotif），只要 knownResults.txt
  local denovo_opt="-nomotif"
  (( ${MOTIF_DENOVO:-0} )) && denovo_opt=""
  (( ${MOTIF_DENOVO:-0} )) || log "只比对已知库（-nomotif）；要 de novo 设 MOTIF_DENOVO=1"

  # -p 不给就按 1 线程跑：homer2 known 单线程扫 1 万 peak 要 7 分钟以上，
  # 三个组分串起来就是半小时。这里显式传 -p，并压到不超过总线程数。
  local mt="${MOTIF_THREADS:-8}"
  (( mt > THREADS )) && mt="$THREADS"
  (( mt < 1 )) && mt=1
  log "HOMER 线程数 -p $mt"

  local frac peaks top outdir
  for frac in $FRACS; do
    peaks="$PEAK/${SAMPLE}.${frac}.filtered.narrowPeak"
    top="$PEAK/${SAMPLE}.${frac}.toppeaks.narrowPeak"
    outdir="$MOTIF/$frac"
    need_file "$peaks" "⑨ 的 ${frac} peak（先跑 --stages peaks）"

    # 不用 sort | head：head 提前退出会让 sort 收到 SIGPIPE，pipefail 下直接中断
    runpipe "sort -k9,9nr '$peaks' -o '$TMPD/.top.$$.sorted' && \
             head -n $MOTIF_TOPN '$TMPD/.top.$$.sorted' > '$top' && \
             rm -f '$TMPD/.top.$$.sorted'"

    runpipe "findMotifsGenome.pl '$top' '$HOMER_GENOME' '$outdir' \
             -size $MOTIF_SIZE -mask -mknown '$HOMER_MOTIFS' $denovo_opt -p $mt \
             > '$LOGS/${SAMPLE}.${frac}.motif.log' 2>&1"
  done

  (( DRY_RUN )) && { mark_done motif; return 0; }

  for frac in $FRACS; do
    local known="$MOTIF/$frac/knownResults.txt"
    local best="" nmotif=""
    if [[ -s "$known" ]]; then
      best=$(awk -F'\t' 'NR==2{print $1}' "$known")
      nmotif=$(awk 'END{print NR-1}' "$known")
    fi

    # 按 fold 排序的整理表，给人看
    [[ -s "$known" ]] && run "$PYTHON" "$PROJECT_ROOT/src/py/motif_top.py" \
        --known "$known" --out "$MOTIF/$frac/top_motifs.tsv" --top 15

    metric "13_motif" "top_motif_${frac}"    "$best"   ""        ""
    metric "13_motif" "motifs_tested_${frac}" "$nmotif" "motifs" ""
  done

  local top1
  top1=$(awk -F'\t' 'NR==2{print $1}' "$MOTIF/all/knownResults.txt" 2>/dev/null || true)
  ok "最强 motif（全部片段）: ${top1:-无}"
  mark_done motif
}
