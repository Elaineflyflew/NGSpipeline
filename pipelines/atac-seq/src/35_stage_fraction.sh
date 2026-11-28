#!/usr/bin/env bash
# =============================================================================
#  src/35_stage_fraction.sh —— ⑦ 按片段长度分档
#
#  为什么要分档：ATAC-seq 的插入片段长度直接反映染色质状态
#    · NFR（< 147 bp）  两个 Tn5 都插在核小体之间的裸 DNA 上 → 开放区信号最强
#    · mono（147-294）  跨一个核小体 → 信号弱，含核小体定位信息
#    · multi（> 294）   跨多个核小体 → 背景
#  好的 ATAC 文库 NFR 应占 40% 以上；NFR 比例低说明细胞核完整性差或过度消化。
#
#  ★ 必须用【已偏移】的 srt2.bam：只有偏移后的 TLEN 才是真实的 Tn5 插入间距。
#    用未偏移的 srt.bam 分档，147 bp 的阈值会错位 9 bp，NFR 档会系统性偏少。
#    分档与 peak calling 必须同源，口径才一致。
# =============================================================================

stage_fraction() {
  section "⑦ 片段长度分档"
  need samtools

  checkfile "$SRT2" "⑥ 产出的 srt2.bam"

  local nfr="$NFR" mono="$MONO" multi="$MULTI"

  # 先统计一遍三档占比（只读一次 BAM），写进 tsv 供报告引用
  if (( ! DRY_RUN )); then
    samtools view -@ "$THREADS" "$SRT2" \
      | awk -v nfrmax="$NFR_MAX" -v monomax="$MONO_MAX" '
          { a = ($9 < 0 ? -$9 : $9)
            if (a == 0)      zero++
            else if (a <  nfrmax)  nfr++
            else if (a <= monomax) mono++
            else                   multi++
            tot++ }
          END {
            printf "total\t%d\n", tot
            printf "nfr\t%d\t%s\n",  nfr,   (tot > 0 ? sprintf("%.1f", 100 * nfr  / tot) : "0")
            printf "mono\t%d\t%s\n", mono,  (tot > 0 ? sprintf("%.1f", 100 * mono / tot) : "0")
            printf "multi\t%d\t%s\n", multi, (tot > 0 ? sprintf("%.1f", 100 * multi / tot) : "0")
            printf "tlen_zero\t%d\n", zero
          }' > "$BAM/${SAMPLE}.fragment_class.tsv"
    log "片段分档: $(awk -F'\t' '$1=="nfr"{print "NFR "$3"%"} $1=="mono"{print "mono "$3"%"} $1=="multi"{print "multi "$3"%"}' "$BAM/${SAMPLE}.fragment_class.tsv" | tr '\n' ' ')"
  fi

  # 三档各写一个 BAM。awk 只判断数据行的 TLEN，header 原样透传
  local tag awkbody
  for tag in nfr mono multi; do
    case "$tag" in
      nfr)   awkbody='{a=($9<0?-$9:$9); if(a<nfrmax && a>0) print}' ;;
      mono)  awkbody='{a=($9<0?-$9:$9); if(a>=nfrmax && a<=monomax) print}' ;;
      multi) awkbody='{a=($9<0?-$9:$9); if(a>monomax) print}' ;;
    esac
    local out
    case "$tag" in
      nfr)   out="$nfr" ;;
      mono)  out="$mono" ;;
      multi) out="$multi" ;;
    esac
    runpipe "samtools view -@ $THREADS -h '$SRT2' \
        | awk -v nfrmax=$NFR_MAX -v monomax=$MONO_MAX '/^@/{print;next} $awkbody' \
        | samtools view -@ $THREADS -b -o '$out' -"
    (( DRY_RUN )) || samtools index -@ "$THREADS" "$out"
  done

  mark_done fraction
}
