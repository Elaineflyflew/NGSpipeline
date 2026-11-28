#!/usr/bin/env bash
# =============================================================================
#  src/53_stage_atacqc.sh —— ⑭ 片段长度分布（内置实现，不依赖 R）
#
#  这一节对应课程里的 ATACseqQC。之所以不用 R 包重写一遍，是因为
#  ATACseqQC 真正不可替代的输出只有两样：
#    · 片段长度分布曲线（看 NFR / 单核小体 / 多核小体的周期峰）
#    · TSS enrichment 分数（⑩ 已经算了）
#  两样都能用 samtools + awk 直接做，省掉一个 R 依赖，也避免版本问题。
#
#  正常的 ATAC 文库在分布图上应该看到：
#    第一个峰 ~50-100 bp（NFR，裸 DNA）
#    第二个峰 ~200 bp（单核小体）
#    第三个峰 ~400 bp（双核小体），之后逐渐衰减
#  只有第一个峰、后面全平 → 过度消化；峰整体右移 → 消化不足。
# =============================================================================

# 用纯 awk 画一张 SVG 折线图，避免引入 gnuplot / R
frag_svg() {
  local tsv="$1" out="$2" nfr="$3" mono="$4" sample="$5"
  awk -F'\t' -v nfr="$nfr" -v mono="$mono" -v sample="$sample" '
    { len[NR] = $1 + 0; cnt[NR] = $2 + 0; if (cnt[NR] > maxc) maxc = cnt[NR]; n = NR; maxl = $1 + 0 }
    END {
      if (n == 0) exit
      W = 900; H = 340; L = 62; R = 18; T = 24; B = 44
      pw = W - L - R; ph = H - T - B

      printf "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 %d %d\" width=\"100%%\" role=\"img\">\n", W, H
      printf "<title>fragment length distribution — %s</title>\n", sample
      printf "<rect width=\"%d\" height=\"%d\" fill=\"#ffffff\"/>\n", W, H

      # y 轴网格 + 刻度
      for (i = 0; i <= 4; i++) {
        yy = T + ph - ph * i / 4
        printf "<line x1=\"%d\" y1=\"%.1f\" x2=\"%d\" y2=\"%.1f\" stroke=\"#E5E7EB\" stroke-width=\"1\"/>\n", L, yy, L + pw, yy
        printf "<text x=\"%d\" y=\"%.1f\" text-anchor=\"end\" font-family=\"system-ui,sans-serif\" font-size=\"11\" fill=\"#6B7280\">%d</text>\n", L - 8, yy + 4, maxc * i / 4
      }

      # x 轴刻度
      for (i = 0; i <= 5; i++) {
        xx = L + pw * i / 5
        printf "<line x1=\"%.1f\" y1=\"%d\" x2=\"%.1f\" y2=\"%d\" stroke=\"#E5E7EB\" stroke-width=\"1\"/>\n", xx, T, xx, T + ph
        printf "<text x=\"%.1f\" y=\"%d\" text-anchor=\"middle\" font-family=\"system-ui,sans-serif\" font-size=\"11\" fill=\"#6B7280\">%d</text>\n", xx, T + ph + 18, maxl * i / 5
      }

      # NFR / mono 分界
      printf "<line x1=\"%.1f\" y1=\"%d\" x2=\"%.1f\" y2=\"%d\" stroke=\"#F59E0B\" stroke-width=\"1\" stroke-dasharray=\"4 3\"/>\n", L + pw * nfr / maxl, T, L + pw * nfr / maxl, T + ph
      printf "<text x=\"%.1f\" y=\"%d\" font-family=\"system-ui,sans-serif\" font-size=\"11\" fill=\"#B45309\">NFR | mono (%d bp)</text>\n", L + pw * nfr / maxl + 4, T + 12, nfr
      printf "<line x1=\"%.1f\" y1=\"%d\" x2=\"%.1f\" y2=\"%d\" stroke=\"#F59E0B\" stroke-width=\"1\" stroke-dasharray=\"4 3\"/>\n", L + pw * mono / maxl, T, L + pw * mono / maxl, T + ph
      printf "<text x=\"%.1f\" y=\"%d\" font-family=\"system-ui,sans-serif\" font-size=\"11\" fill=\"#B45309\">mono | multi (%d bp)</text>\n", L + pw * mono / maxl + 4, T + 28, mono

      # 折线
      printf "<polyline fill=\"none\" stroke=\"#2563EB\" stroke-width=\"1.8\" points=\""
      for (i = 1; i <= n; i++) {
        printf "%.1f,%.1f ", L + pw * len[i] / maxl, T + ph - ph * cnt[i] / maxc
      }
      printf "\"/>\n"

      # 轴与标题
      printf "<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"#9CA3AF\" stroke-width=\"1\"/>\n", L, T + ph, L + pw, T + ph
      printf "<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"#9CA3AF\" stroke-width=\"1\"/>\n", L, T, L, T + ph
      printf "<text x=\"%d\" y=\"%d\" font-family=\"system-ui,sans-serif\" font-size=\"13\" font-weight=\"500\" fill=\"#111827\">片段长度分布 · %s</text>\n", L, 16, sample
      printf "<text x=\"%d\" y=\"%d\" text-anchor=\"middle\" font-family=\"system-ui,sans-serif\" font-size=\"12\" fill=\"#6B7280\">片段长度 (bp)</text>\n", L + pw / 2, H - 8
      printf "</svg>\n"
    }' "$tsv" > "$out"
}

stage_atacqc() {
  if (( ! RUN_ATACQC )); then
    log "⑭ 跳过片段分布（RUN_ATACQC=0）"
    mark_done atacqc
    return 0
  fi

  section "⑭ 片段长度分布"
  need samtools
  checkfile "$SRT2" "⑥ 产出的 srt2.bam"

  local tsv="$AQC/${SAMPLE}.fragment_length.tsv"

  if (( ! DRY_RUN )); then
    # 5 bp 一个 bin，只统计到 1000 bp（再长基本是噪声）
    samtools view -@ "$THREADS" "$SRT2" \
      | awk '{ a = ($9 < 0 ? -$9 : $9); if (a > 0 && a <= 1000) h[int(a / 5) * 5]++ }
             END { for (k = 0; k <= 1000; k += 5) if (k in h) printf "%d\t%d\n", k, h[k] }' \
      > "$tsv"

    frag_svg "$tsv" "$AQC/${SAMPLE}.fragment_length.svg" "$NFR_MAX" "$MONO_MAX" "$SAMPLE"

    # 主峰位置：判断"像不像 ATAC"的一个数
    local mode
    mode=$(sort -k2,2nr "$tsv" | head -1 | cut -f1)
    printf 'fragment_mode_bp\t%s\n' "$mode" >> "$AQC/${SAMPLE}.fragment_class_extra.tsv"
    log "片段长度主峰: ${mode} bp（预期 50-100 bp）"
  fi

  mark_done atacqc
}
