#!/usr/bin/env bash
################################################################################
# 离线单测：不依赖任何生信工具，只用合成数据验证 src/py/ 下的解析与计算逻辑。
#
#   bash tests/run_tests.sh
#
# 覆盖：
#   · fastp_stats    JSON 字段抽取
#   · tss_prep       GTF → TSS（正负链取点方向）
#   · tss_score      computeMatrix 矩阵 → TSS 富集分数
#   · annot_dist     HOMER annStats → 分布表
#   · motif_top      knownResults → top motif（fold 排序）
#   · collect_metrics  阈值判读 PASS/WARN
#   · make_report    HTML 生成
#
# 不在覆盖范围：frag_bin.py（需要真 BAM + samtools）、atacqc.R（需要 Bioconductor）。
################################################################################

set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$HERE")"
SCRIPTS="$ROOT/src/py"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

PASS=0
FAIL=0

check() {  # check <描述> <期望出现的字符串> <文件>
    local desc="$1" want="$2" file="$3"
    if [[ ! -s "$file" ]]; then
        printf '  \033[31mFAIL\033[0m %s —— 输出文件为空\n' "$desc"
        FAIL=$((FAIL + 1)); return
    fi
    if grep -q -- "$want" "$file"; then
        printf '  \033[32mPASS\033[0m %s\n' "$desc"
        PASS=$((PASS + 1))
    else
        printf '  \033[31mFAIL\033[0m %s —— 没找到 "%s"\n' "$desc" "$want"
        printf '        实际内容：\n'; sed 's/^/          /' "$file" | head -8
        FAIL=$((FAIL + 1))
    fi
}

echo "== 准备合成数据 =="

cat > "$TMP/fastp.json" <<'EOF'
{
  "summary": {
    "before_filtering": {"total_reads": 40000000, "q30_rate": 0.931, "gc_content": 0.412},
    "after_filtering":  {"total_reads": 39600000, "q30_rate": 0.962}
  },
  "filtering_result": {"low_quality_reads": 120000, "too_short_reads": 280000},
  "duplication": {"rate": 0.121},
  "adapter_cutting": {"adapter_trimmed_reads": 1500000}
}
EOF

cat > "$TMP/test.gtf" <<'EOF'
chr1	HAVANA	gene	1000	2000	.	+	.	gene_id "ENSG1"; gene_name "A"; gene_type "protein_coding";
chr1	HAVANA	gene	5000	4000	.	-	.	gene_id "ENSG2"; gene_name "B"; gene_type "protein_coding";
chr1	HAVANA	gene	9000	9500	.	+	.	gene_id "ENSG3"; gene_name "C"; gene_type "lncRNA";
EOF

# 真实形状：bin size 25、±2000 bp → 160 个 bin。中心 3 个 bin 信号 5，其余 1
# → TSS 分数应为 5.0
python3 - "$TMP/matrix.gz" <<'PY'
import gzip, json, sys
n = 160
hdr = {"sample_labels": ["x.bw"], "sample_boundaries": [3, 3 + n],
       "bin size": 25, "upstream": 2000, "downstream": 2000}
with gzip.open(sys.argv[1], "wt") as fh:
    fh.write("@" + json.dumps(hdr) + "\n")
    for r in range(10):
        vals = [5 if i in (79, 80, 81) else 1 for i in range(n)]
        fh.write("chr1\t%d\t%d\t%s\n" % (r * 100, r * 100 + 1,
                 "\t".join(str(v) for v in vals)))
PY

cat > "$TMP/annostats.txt" <<'EOF'
Annotation	Number of Peaks	% of total
Promoter (<=1kb)	20732	27.96
Intron (1-5kb)	12000	16.19
Distal Intergenic	25000	33.72
Total	74134	100.00
EOF

cat > "$TMP/knownResults.txt" <<'EOF'
Motif Name	Consensus	P-value	Log P-value	q-value (Benjamini)	# of Target Sequences with Motif(of 10000)	% of Target Sequences with Motif	# of Background Sequences with Motif(of 45678)	% of Background Sequences with Motif
NFY(CCAAT)/HepG2	AGCCAATCAG	1e-300	-690.8	0	3000	30.00%	500	1.10%
CTCF(Zf)/CD4+-CTCF-ChIP-Seq	CCGCACAGGGGGCGCT	1e-500	-1151.0	0	4800	48.00%	300	0.66%
EOF

mkdir -p "$TMP/metrics"
printf '02_fastp\treads_retained\t0.99\tratio\t>0.95\n'  > "$TMP/metrics/02_fastp.tsv"
printf '03_bowtie2\tchrm_fraction\t0.2734\tratio\t<0.05\n' > "$TMP/metrics/03_bowtie2.tsv"
printf '07_fragbin\tnfr_fraction\t0.41\tratio\t>0.30\n'   > "$TMP/metrics/07_fragbin.tsv"
printf '10_frip\tfrip_all\t0.5388\tratio\t>0.20\n'       > "$TMP/metrics/10_frip.tsv"
printf '11_coverage\ttss_enrichment_score\t8.20\tfold\t>6\n' > "$TMP/metrics/11_coverage.tsv"

echo
echo "== fastp_stats =="
python3 "$SCRIPTS/fastp_stats.py" --json "$TMP/fastp.json" --out "$TMP/fastp.tsv" >/dev/null
check "reads 对数 = 2000 万" "read_pairs_before_trim	20000000" "$TMP/fastp.tsv"
check "保留率 0.99"          "reads_retained	0.99"             "$TMP/fastp.tsv"

echo
echo "== tss_prep =="
python3 "$SCRIPTS/tss_prep.py" --gtf "$TMP/test.gtf" --out "$TMP/tss.bed" >/dev/null
check "正链 TSS = start-1" "chr1	999	1000	A	0	+"   "$TMP/tss.bed"
check "负链 TSS = end"     "chr1	4000	4001	B	0	-"  "$TMP/tss.bed"
if [[ "$(wc -l < "$TMP/tss.bed")" -eq 2 ]]; then
    printf '  \033[32mPASS\033[0m lncRNA 被过滤（只剩 2 个 TSS）\n'; PASS=$((PASS + 1))
else
    printf '  \033[31mFAIL\033[0m lncRNA 没被过滤\n'; FAIL=$((FAIL + 1))
fi

echo
echo "== tss_score =="
python3 "$SCRIPTS/tss_score.py" --matrix "$TMP/matrix.gz" --flank 2000 --window 50 \
    --out "$TMP/tss_score.tsv" >/dev/null
check "TSS 富集分数 = 5.000" "tss_enrichment_score	5.000" "$TMP/tss_score.tsv"

echo
echo "== annot_dist =="
python3 "$SCRIPTS/annot_dist.py" --annstats "$TMP/annostats.txt" --out "$TMP/annot.tsv" >/dev/null
check "解析出 Distal Intergenic" "Distal Intergenic	25000	33.72" "$TMP/annot.tsv"

echo
echo "== motif_top =="
python3 "$SCRIPTS/motif_top.py" --known "$TMP/knownResults.txt" --out "$TMP/motif.tsv" >/dev/null
check "按 fold 排序，CTCF 第一" "1	CTCF" "$TMP/motif.tsv"
check "fold 计算正确（48/0.66≈72.7）" "72.7" "$TMP/motif.tsv"

echo
echo "== collect_metrics =="
python3 "$SCRIPTS/collect_metrics.py" --metrics "$TMP/metrics" --out "$TMP/summary.tsv" >/dev/null
check "FRiP 0.54 判 PASS"   "0.5388	ratio	>0.20	PASS"  "$TMP/summary.tsv"
check "chrM 27% 判 WARN"    "0.2734	ratio	<0.05	WARN"  "$TMP/summary.tsv"

echo
echo "== make_report =="
mkdir -p "$TMP/report_dir"
cp "$TMP/summary.tsv" "$TMP/report_dir/summary.tsv"
python3 "$SCRIPTS/make_report.py" --summary "$TMP/report_dir/summary.tsv" \
    --sample TEST --outdir "$TMP/report_dir" >/dev/null
check "生成了 HTML 报告" "<!DOCTYPE html>" "$TMP/report_dir/report.html"
check "报告里带出了结论"  "线粒体占比"      "$TMP/report_dir/report.html"
check "三项齐全时给合看结论" "三指标合看"   "$TMP/report_dir/report.html"
check "TSS 8.20× 判达标"     "TSS 富集 8.20× 达标" "$TMP/report_dir/report.html"

echo
printf '== 结果：%d 通过 / %d 失败 ==\n' "$PASS" "$FAIL"
[[ $FAIL -eq 0 ]] || exit 1
