#!/usr/bin/env bash
################################################################################
# tests/smk_dryrun.sh —— Snakemake 工作流的离线 DAG 校验
#
#   bash tests/smk_dryrun.sh
#
# 干什么：造一套假数据（空 FASTQ、空参考、两个 HOMER 假命令），让
# snakemake 只构建 DAG 不执行，确认：
#   · 所有 rule 的 input/output/log/benchmark 通配符自洽
#   · 没有 AmbiguousRuleException（同一个文件能被两条 rule 生产）
#   · 没有残留的 {res(...)} / {SCRIPTS} 之类非 Snakemake 占位符
#   · 15 步对应的 rule 全部进了 DAG
#
# 为什么值得单独写：这套流程真正的执行需要 bowtie2/samtools/MACS2/HOMER，
# 本机装不了；但 DAG 能不能建起来跟工具有没有装无关，能在提交前就查出来。
#
# 没装 snakemake 时直接跳过（不算失败）。
################################################################################

set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$HERE")"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

PASS=0
FAIL=0
ok()   { printf '  \033[32mPASS\033[0m %s\n' "$1"; PASS=$((PASS + 1)); }
bad()  { printf '  \033[31mFAIL\033[0m %s\n' "$1"; FAIL=$((FAIL + 1)); }

if ! command -v snakemake >/dev/null 2>&1; then
    printf '\033[33m[SKIP]\033[0m 没有 snakemake，跳过 DAG 校验\n'
    printf '       装：mamba install -c bioconda -c conda-forge snakemake-minimal\n'
    exit 0
fi

echo "== 造假数据 =="
mkdir -p "$TMP"/{bin,rawdata,refs/bt2,out}

# 空 gzip：只为让路径存在，snakemake 不做内容校验
for i in 1 2; do : > "$TMP/rawdata/demo_${i}.fastq.gz"; done

: > "$TMP/refs/genes.gtf"
: > "$TMP/refs/hg38.fa"
: > "$TMP/refs/blacklist.bed"
: > "$TMP/refs/known.motifs"
for ext in 1.bt2 2.bt2 3.bt2 4.bt2 rev.1.bt2 rev.2.bt2; do
    : > "$TMP/refs/bt2/hg38.$ext"
done

# HOMER 假命令：只为了让 homer_available() 返回 True，把 ⑫⑬ 也拉进 DAG
for t in annotatePeaks.pl findMotifsGenome.pl; do
    printf '#!/bin/sh\nexit 0\n' > "$TMP/bin/$t"
    chmod +x "$TMP/bin/$t"
done

cat > "$TMP/organism.yaml" <<EOF
genome_size: 2900338458
bowtie2_index: "$TMP/refs/bt2/hg38"
genome_fasta: "$TMP/refs/hg38.fa"
genes_gtf: "$TMP/refs/genes.gtf"
blacklist_bed: "$TMP/refs/blacklist.bed"
homer_genome: "hg38"
homer_known_motifs: "$TMP/refs/known.motifs"
ignoreForNormalization: ""
EOF

# 只写需要覆盖的键：config/config.yaml 仍然生效，这个文件叠在上面
cat > "$TMP/config.yaml" <<EOF
organism_file: "$TMP/organism.yaml"
threads: 2
samples:
  sub20M.SRR13849313:
    r1: "$TMP/rawdata/demo_1.fastq.gz"
    r2: "$TMP/rawdata/demo_2.fastq.gz"
EOF

echo
echo "== 构建 DAG =="
LOG="$TMP/dag.log"
PATH="$TMP/bin:$PATH" snakemake \
    -s "$ROOT/workflow/Snakefile" \
    --configfile "$TMP/config.yaml" \
    -n -j 2 --nolock --printshellcmds \
    > "$LOG" 2>&1
RC=$?

if [[ $RC -eq 0 ]]; then
    ok "snakemake -n 正常退出"
else
    bad "snakemake -n 退出码 $RC"
    echo "        最后 25 行："
    tail -25 "$LOG" | sed 's/^/        /'
fi

# ---------------------------------------------------------------------------
# 逐项断言
# ---------------------------------------------------------------------------
if grep -q "AmbiguousRuleException" "$LOG"; then
    bad "存在规则歧义（AmbiguousRuleException）"
    grep -A3 "AmbiguousRuleException" "$LOG" | sed 's/^/        /'
else
    ok "没有规则歧义"
fi

# 非 Snakemake 占位符漏进 shell 的典型症状
STALE='\{(res\(|SCRIPTS\}|RES\}|LOGS\}|BENCH\}|SAMPLES\[|PROJECT_ROOT\})'
if grep -qE "$STALE" "$LOG"; then
    bad "shell 命令里有未展开的 Python 表达式"
    grep -nE "$STALE" "$LOG" | head -5 | sed 's/^/        /'
else
    ok "没有残留的 Python 占位符"
fi

# 15 步对应的 rule 是否都进了 DAG
expect_rules=(
    fastqc_raw fastp fastqc_clean
    bowtie2 bam_filter atac_shift resort frag_bin
    macs2 blacklist
    frip bam_coverage tss_prep tss_matrix tss_profile tss_score
    homer_annotate annot_dist
    homer_motif motif_top
    atacqc igv_snapshot
    collect_metrics report
)
missing=()
for r in "${expect_rules[@]}"; do
    grep -qE "^\s*$r\s+[0-9]+$" "$LOG" || missing+=("$r")
done
if [[ ${#missing[@]} -eq 0 ]]; then
    ok "15 步对应 ${#expect_rules[@]} 条 rule 全部在 DAG 里"
else
    bad "DAG 里缺 rule：${missing[*]}"
fi

# 三条 peak calling 路线（all / nfr / mono）都跑了
for frac in all nfr mono; do
    if grep -q "${frac}.filtered.narrowPeak" "$LOG"; then
        ok "路线 $frac 有 peak"
    else
        bad "路线 $frac 缺 peak"
    fi
done

# 三条路线的 FRiP 都算了
for frac in all nfr mono; do
    if grep -q "/${frac}.frip.tsv" "$LOG"; then
        ok "路线 $frac 有 FRiP"
    else
        bad "路线 $frac 缺 FRiP"
    fi
done

# 关键参数是否真的出现在命令行里（防止 config 里改了但没接上）
grep -q -- "--ATACshift" "$LOG"     && ok "⑤ alignmentSieve 带 --ATACshift" || bad "⑤ 缺 --ATACshift"
grep -q -- "-f BAMPE" "$LOG"        && ok "⑧ MACS2 用 BAMPE 模式"        || bad "⑧ 不是 BAMPE 模式"
grep -q -- "--nomodel" "$LOG"       && ok "⑧ MACS2 带 --nomodel"         || bad "⑧ 缺 --nomodel"
grep -q -- "intersect -v" "$LOG"    && ok "⑨ 黑名单用 intersect -v"      || bad "⑨ 黑名单用法不对"
grep -q -- "--normalizeUsing CPM" "$LOG" && ok "⑪ bamCoverage CPM 归一"  || bad "⑪ 缺 CPM 归一"

# MACS2 绝对不能带 --shift/--extsize（BAMPE 模式下会被静默忽略，见 README 决策 3）
if grep -q -- "--shift\|--extsize" "$LOG"; then
    bad "⑧ 出现了 --shift/--extsize（BAMPE 模式下静默无效）"
else
    ok "⑧ 没有误传 --shift/--extsize"
fi

echo
printf '== 结果：%d 通过 / %d 失败 ==\n' "$PASS" "$FAIL"
[[ $FAIL -eq 0 ]] || exit 1
