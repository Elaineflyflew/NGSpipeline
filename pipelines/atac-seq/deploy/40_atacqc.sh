#!/usr/bin/env bash
# =============================================================================
#  deploy/40_atacqc.sh —— ATACseqQC（Bioconductor）
#
#  这一步是可选的：⑭ 的作用是"换个实现再算一遍 TSSE"，和 ⑪ 交叉验证。
#  装不上不影响主链——⑦⑪ 已经给出同口径的片段分布和 TSS 富集数值。
#
#  依赖链比较长（GenomicAlignments / Rsamtools / TxDb.Hsapiens.UCSC.hg38.knownGene），
#  首次编译要 20–40 分钟，别以为卡死了。
# =============================================================================

deploy_atacqc() {
  section "ATACseqQC（R / Bioconductor）"
  need Rscript

  if (( ! DRY_RUN )) && Rscript -e \
      'if (!requireNamespace("ATACseqQC", quietly=TRUE)) quit(status=1)' >/dev/null 2>&1; then
    log "ATACseqQC 已安装，跳过"
    return 0
  fi

  if (( DRY_RUN )); then
    printf '[dry-run] Rscript -e "install.packages(\"BiocManager\"); BiocManager::install(\"ATACseqQC\")"\n'
    return 0
  fi

  log "安装 BiocManager + ATACseqQC（首次约 20–40 分钟）"
  Rscript -e '
    if (!requireNamespace("BiocManager", quietly=TRUE))
      install.packages("BiocManager", repos="https://cloud.r-project.org")
    BiocManager::install(c("ATACseqQC", "TxDb.Hsapiens.UCSC.hg38.knownGene"),
                         ask=FALSE, update=FALSE)
  ' || die "ATACseqQC 安装失败（可跳过：⑭ 不是主链必需）"

  ok "ATACseqQC 就绪"
}
