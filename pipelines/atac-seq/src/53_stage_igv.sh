#!/usr/bin/env bash
# =============================================================================
#  src/53_stage_igv.sh —— ⑮ 位点截图
#
#  用 batch 脚本而不是手动点：可复现、可进 CI、可写进报告附录。
#  服务器没有图形环境时用 xvfb-run；没有 IGV 就退回 pyGenomeTracks 出同一区域。
#
#  读图顺序（课件里那张 nf-core 官方截图的读法）：
#    先看 BED 条（peak 落在哪）→ 再看 bigWig 形态（有没有富集）
#    → 最后看底部基因模型（落在启动子还是远端增强子）。
#  比较幅度之前必须先确认归一化方式——这里统一是 CPM。
# =============================================================================

stage_igv() {
  section "⑮ 位点截图（${IGV_LOCUS}）"
  need_file "$COV/${SAMPLE}.bamCoverage.bw" "⑪ 的 bigWig（先跑 --stages coverage）"

  local peaks="$PEAK/${SAMPLE}.all.filtered.narrowPeak"
  local png="$IGV/${SAMPLE}.MYC.png"
  local batch="$IGV/${SAMPLE}.igv_batch.txt"
  local tracks="$IGV/${SAMPLE}.tracks.ini"

  # batch 文件无论如何都生成：没有图形环境时用户可拷到本地跑
  # genome 指向本地 fasta 而不是写 hg38：写 hg38 时 IGV 会去它的服务器下载
  # 2.6 G 的基因组（离线/内网环境直接卡死），本地 fasta 有 .fai 就能直接用。
  local genome_arg="hg38"
  [[ -n "$GENOME" && -s "$GENOME" ]] && genome_arg="$GENOME"

  # 基因模型：不加载的话截图里只有信号和 peak 条，看不出 peak 落在启动子还是
  # 远端增强子（而这正是"读图三步"的最后一步）。
  # 但**不能直接 load 未建索引的 GTF**：gencode v43 有 150 万行，IGV 会把整个
  # 文件读进内存，实测 10 分钟还没画完（timeout 600 直接把它砍了）。
  # 所以 auto 模式只认"已排序 + tabix 建过索引"的 .gtf.gz —— 那样 IGV 只读目标
  # 区间，秒级出图。没索引就不加载，日志里说明怎么办。
  #   IGV_GENE_TRACK=auto（默认）| 0（关）| 其他值（当成路径）
  local ann_dir="$PROJECT_ROOT/reference/annotation"
  [[ -n "${GTF:-}" ]] && ann_dir=$(dirname "$GTF")
  local gene_track=""
  case "${IGV_GENE_TRACK:-auto}" in
    0|off|no) gene_track="" ;;
    auto)
      local g
      for g in "$ann_dir"/*.gtf.gz "$ann_dir"/*.gff3.gz; do
        [[ -s "$g" && -s "$g.tbi" ]] && { gene_track="$g"; break; }
      done
      if [[ -z "$gene_track" ]]; then
        log "没找到已建索引的 *.gtf.gz，截图不带基因模型（做法见 RUNBOOK 第 7 节）"
      fi
      ;;
    *) gene_track="$IGV_GENE_TRACK" ;;
  esac
  local load_gene=""
  [[ -n "$gene_track" ]] && load_gene="load $gene_track"
  [[ -n "$gene_track" ]] && log "基因模型: $(basename "$gene_track")"

  runpipe "cat > '$batch' <<'EOF'
new
genome $genome_arg
load $COV/${SAMPLE}.bamCoverage.bw
load $peaks
$load_gene
goto $IGV_LOCUS
sort position
snapshot $png
exit
EOF"

  local igv_bin="" pgt=""
  igv_bin=$(find_tool igv 2>/dev/null || true)
  pgt=$(find_tool pyGenomeTracks 2>/dev/null || true)

  # IGV 的 batch 模式偶尔会卡在"加载基因组/渲染"上不退出（尤其是 xvfb 下），
  # 所以套一层 timeout：宁可没图，也别让整条链挂在这里。
  local to=""
  if find_tool timeout >/dev/null 2>&1; then
    to="timeout ${IGV_TIMEOUT:-600}"
  fi

  if [[ -n "$igv_bin" ]]; then
    # 用 try_runpipe：IGV 起不来时退出码可能是 1，不该把整条链带崩
    if find_tool xvfb-run >/dev/null 2>&1; then
      try_runpipe "$to xvfb-run -a '$igv_bin' -b '$batch' > '$LOGS/${SAMPLE}.igv.log' 2>&1"
    else
      try_runpipe "$to '$igv_bin' -b '$batch' > '$LOGS/${SAMPLE}.igv.log' 2>&1"
    fi
  elif [[ -n "$pgt" ]]; then
    log "没有 IGV，改用 pyGenomeTracks"
    runpipe "cat > '$tracks' <<'EOF'
[bigwig]
file = $COV/${SAMPLE}.bamCoverage.bw
title = $SAMPLE (CPM)
height = 6
color = #8B0000

[peaks]
file = $peaks
title = MACS2 peaks
height = 2
color = #1f4e79
EOF"
    try_runpipe "$to '$pgt' --tracks '$tracks' --region '$IGV_LOCUS' -o '$png' \
             >> '$LOGS/${SAMPLE}.igv.log' 2>&1"
  else
    if (( DRY_RUN )); then
      printf '[dry-run] 既没有 igv 也没有 pyGenomeTracks，只会生成 batch 文件\n'
      mark_done igv
      return 0
    fi
    warn "既没有 igv 也没有 pyGenomeTracks；batch 已生成，可拷到本地跑"
    metric "15_igv" "snapshot_generated" "0" "" ""
    mark_done igv
    return 0
  fi

  (( DRY_RUN )) && { mark_done igv; return 0; }

  # IGV 退出码为 0 不代表真的截了图（窗口没起来就静默退出了），
  # 所以以"PNG 是否真的存在"为准。
  if [[ -s "$png" ]]; then
    metric "15_igv" "snapshot_generated" "1" "" ""
    ok "截图: $png"
  else
    metric "15_igv" "snapshot_generated" "0" "" ""
    warn "IGV 没产出 PNG（看 $LOGS/${SAMPLE}.igv.log）；batch 文件已生成，可拷到本地跑"
  fi
  mark_done igv
}
