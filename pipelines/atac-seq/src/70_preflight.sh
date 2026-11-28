#!/usr/bin/env bash
# =============================================================================
#  src/70_preflight.sh —— 跑数据前的体检
#
#  分两半：
#    工具   只检查"当前 STAGES 里真正会用到的"，跑 --stages peaks 时不该因为
#           没装 ATACseqQC 就拦下来。
#    参考   基因组 / bowtie2 索引 / GTF / 黑名单，缺一个后面都跑不通。
#
#  --check 时缺了只 WARN 并计数，最后统一汇报；正常运行时第一个缺就退出。
# =============================================================================

# 每个阶段真正会用到的命令。按阶段列而不是给一条"核心命令"清单，
# 是因为 --stages qc 时不该因为没装 bowtie2 就拦下来。
stage_tools() {
  case "$1" in
    qc)       printf 'fastqc' ;;
    trim)     printf 'fastp fastqc' ;;
    align)    printf 'bowtie2 samtools' ;;
    bam)      printf 'samtools' ;;
    frag)     printf 'samtools python3' ;;
    peaks)    printf 'macs2 samtools bedtools' ;;
    frip)     printf 'samtools bedtools python3' ;;
    coverage) printf 'samtools bamCoverage computeMatrix plotProfile plotFingerprint python3' ;;
    annot)    printf 'annotatePeaks.pl python3' ;;
    motif)    printf 'findMotifsGenome.pl python3' ;;
    atacqc)   printf 'Rscript samtools' ;;
    igv)      printf 'samtools' ;;
    report)   printf 'python3 multiqc' ;;
    *)        printf '' ;;
  esac
}

# 参考数据缺失的三种处理：dry-run 只提示 / --check 记账 / 正常运行退出
ref_missing() {
  local what="$1" hint="$2"
  if (( DRY_RUN )); then
    printf '[dry-run] 缺 %s（实际运行会报错）\n' "$what"
    SOFT_MISSING=$(( SOFT_MISSING + 1 ))
    return 0
  fi
  if (( CHECK_ONLY )); then
    warn "缺 $what"
    MISSING_TOOLS=$(( MISSING_TOOLS + 1 ))
    return 0
  fi
  die "缺 $what —— $hint"
}

preflight() {
  section "体检：工具与参考数据"

  # 按阶段取并集（13 个阶段有大量重复命令，用 seen 去重，
  # 否则缺一个 samtools 会报六遍）
  local t s seen=" "
  for s in ${STAGES//,/ }; do
    for t in $(stage_tools "$s"); do
      case "$seen" in *" $t "*) continue ;; esac
      seen="$seen$t "

      # multiqc 是可选聚合，缺了 report 仍能出（自己生成的 HTML 才是主交付物）
      if [[ "$t" == "multiqc" ]] && ! find_tool multiqc >/dev/null 2>&1; then
        warn "没有 multiqc，报告里的聚合页会跳过（不影响 report.html）"
        continue
      fi
      need "$t"
    done
  done

  # ---- 参考数据 ----
  # dry-run 只提示、--check 记一笔、正常运行直接退出
  #
  # 只有 align 及之后的阶段才读基因组，只跑 qc/trim 时不该被"没有参考"拦住。
  # 这样在没装参考的机器上也能先验一遍 ①② 两步。
  local needs_genome=0 needs_bt2=0
  for s in ${STAGES//,/ }; do
    case "$s" in
      align) needs_genome=1; needs_bt2=1 ;;
      bam|frag|peaks|frip|coverage|annot|motif|atacqc|igv) needs_genome=1 ;;
    esac
  done

  if (( ! needs_genome )); then
    log "本次只跑 ${STAGES}，不读基因组"
  elif [[ -z "$GENOME" || ! -s "$GENOME" ]]; then
    ref_missing "基因组 fasta（reference/genome/）" "跑 deploy/deploy.sh genome，或用 -g 指定"
  else
    log "基因组:   $GENOME"
  fi

  if (( needs_bt2 )); then
    if [[ -z "$BT2" ]]; then
      ref_missing "bowtie2 索引（reference/genome/bt2/）" "跑 deploy/deploy.sh bt2"
    else
      # 只查 .1.bt2 不够：缺 rev.1/rev.2 时 bowtie2 会 segfault 而不是报错
      need_bt2 "$BT2"
      log "bt2 索引: $BT2"
    fi
  fi

  # GTF 只有 ⑪ 的 TSS 富集要用；黑名单只有 ⑧⑨ 和 ⑪ 要用。用不上就别念
  if stage_enabled coverage; then
    if [[ -z "$GTF" ]]; then
      warn "没找到 GTF，⑪ 的 TSS 富集会被跳过（deploy/deploy.sh annotation）"
    else
      log "GTF:      $GTF"
    fi
  fi

  if stage_enabled peaks || stage_enabled coverage; then
    if [[ -z "$BLACKLIST" ]]; then
      warn "没找到黑名单 BED，⑨ 会退化成“原样拷贝”（deploy/deploy.sh blacklist）"
    else
      log "黑名单:   $BLACKLIST"
    fi
  fi

  # ---- 两条可选链路 ----
  if stage_enabled igv; then
    if find_tool igv >/dev/null 2>&1 || find_tool pyGenomeTracks >/dev/null 2>&1; then
      log "IGV 路线: $(find_tool igv 2>/dev/null || find_tool pyGenomeTracks)"
    else
      warn "igv 和 pyGenomeTracks 都没有，⑮ 只会生成 batch 文件"
    fi
  fi

  if (( MISSING_TOOLS > 0 )); then
    warn "共缺 ${MISSING_TOOLS} 项（见上）"
  elif (( SOFT_MISSING > 0 )); then
    warn "dry-run：另有 ${SOFT_MISSING} 项参考数据缺失，实际运行会在这里停下（见上）"
  else
    ok "工具与参考数据齐全"
  fi
}
