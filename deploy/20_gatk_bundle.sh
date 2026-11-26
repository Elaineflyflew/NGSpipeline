#!/usr/bin/env bash
# =============================================================================
#  deploy/20_gatk_bundle.sh —— GATK 已知位点资源（hg38）
#
#  四个文件各有用途，缺一个就少一项能力：
#    dbsnp_146.hg38               BQSR 必需，也是变异注释里的 rsID 来源
#    Mills_and_1000G_gold_standard BQSR 的 indel 已知位点
#    1000G_phase1.snps.high_conf   BQSR 的 SNP 已知位点
#    hapmap_3.3.hg38              FilterVariantTranches 分层必需
#
#  全部约 1.5G，比参考基因组小得多。
#
#  已经手动准备好 bundle 的话（解压的 .vcf + .idx 也算），本脚本不会重复下载：
#  只要同名 .vcf 或 .vcf.gz 存在就跳过。
# =============================================================================
set -euo pipefail

BROAD_BASE="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0"

deploy_bundle() {
  mkdir -p "$BUNDLE_DIR"
  check_space 3 "$BUNDLE_DIR"

  local n
  for n in $GATK_BUNDLE_FILES; do
    if have_bundle_vcf "$n"; then
      log "已存在，跳过: ${n}（.vcf 或 .vcf.gz 任一即可）"
      continue
    fi
    dl "$BROAD_BASE/$n.vcf.gz"     "$BUNDLE_DIR/$n.vcf.gz"
    dl "$BROAD_BASE/$n.vcf.gz.tbi" "$BUNDLE_DIR/$n.vcf.gz.tbi"
  done

  if (( DRY_RUN )); then
    log "（dry-run）bundle 将放在 $BUNDLE_DIR"
    return 0
  fi

  # 下载完确认能被解析：文件截断时这里会立刻暴露，而不是等 BQSR 跑到一半
  local f ok=0 bad=0
  if command -v bcftools >/dev/null 2>&1; then
    for n in $GATK_BUNDLE_FILES; do
      f=$(bundle_file "$n") || continue
      if bcftools view -h "$f" >/dev/null 2>&1; then
        log "OK  $(basename "$f")"
        ok=$(( ok + 1 ))
      else
        warn "文件可能损坏: $(basename "$f")（删掉该文件重跑本脚本可重下）"
        bad=$(( bad + 1 ))
      fi
    done
  fi

  if (( bad > 0 )); then
    warn "$BUNDLE_DIR 里有 ${bad} 个文件无法解析"
  fi
  log "bundle 就绪: ${BUNDLE_DIR}（${ok}/4 个已知位点）"
}

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
  source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
  deploy_bundle "$@"
fi
