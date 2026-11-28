#!/usr/bin/env bash
# =============================================================================
#  deploy/lib.sh —— deploy 的公共层
#
#  和主程序共用 src/00_common.sh 与 config/pipeline.conf，
#  这样 doctor 看到的工具清单和 pipeline 运行时用的完全一致，
#  不会出现"doctor 说 OK、跑起来却找不到"这种事。
# =============================================================================

PROJECT_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
export PROJECT_ROOT

# deploy 是独立入口，没有主程序的默认值，这里补上
: "${DRY_RUN:=0}"
: "${CHECK_ONLY:=0}"
: "${FORCE:=0}"
: "${STAGES:=}"
: "${STATUS:=/tmp}"
: "${TMPD:=/tmp}"

# shellcheck source=../src/00_common.sh
source "$PROJECT_ROOT/src/00_common.sh"
# shellcheck source=../config/pipeline.conf
source "$PROJECT_ROOT/config/pipeline.conf"
# shellcheck source=../src/70_preflight.sh
source "$PROJECT_ROOT/src/70_preflight.sh"

# conda 包名和可执行文件名不一定相同
pkg_for() {
  case "$1" in
    macs2)     printf 'macs2' ;;
    bedtools)  printf 'bedtools' ;;
    bowtie2)   printf 'bowtie2' ;;
    fastqc)    printf 'fastqc' ;;
    fastp)     printf 'fastp' ;;
    samtools)  printf 'samtools' ;;
    alignmentSieve|bamCoverage|computeMatrix|plotProfile|plotHeatmap|plotEnrichment|plotFingerprint)
               printf 'deeptools' ;;
    annotatePeaks.pl|findMotifsGenome.pl) printf 'homer' ;;
    *)         printf '%s' "$1" ;;
  esac
}

# 哪些工具在 conda 里有、哪些必须手动装
conda_has() {
  case "$1" in
    fastqc|fastp|bowtie2|samtools|macs2|bedtools|multiqc|deeptools|homer) return 0 ;;
    *) return 1 ;;
  esac
}
