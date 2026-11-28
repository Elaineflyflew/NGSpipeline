#!/usr/bin/env bash
# =============================================================================
#  tests/smoke_test.sh —— 冒烟测试（不需要任何生信工具）
#
#  只验证"脚手架"是否正确：参数解析、路径派生、阶段跳过、dry-run 不落标记、
#  ​--stages 真的只跑指定阶段。真实数据必须上服务器跑，这里连 bowtie2 都不需要装。
#
#    bash tests/smoke_test.sh
# =============================================================================
set -uo pipefail

PROJECT_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
PIPE="$PROJECT_ROOT/bin/atac"
DEPLOY="$PROJECT_ROOT/deploy/deploy.sh"

T=$(mktemp -d "${TMPDIR:-/tmp}/atac-smoke.XXXXXX")
trap 'rm -rf "$T"' EXIT

PASS=0; FAIL=0; FAILED_NAMES=()

ok()   { PASS=$((PASS+1)); printf '  \033[32mPASS\033[0m %s\n' "$1"; }
bad()  { FAIL=$((FAIL+1)); FAILED_NAMES+=("$1"); printf '  \033[31mFAIL\033[0m %s\n' "$1"; [[ -n "${2:-}" ]] && printf '       %s\n' "$2"; }

# expect <名称> <期望退出码> <命令...>
expect() {
  local name="$1" want="$2"; shift 2
  local out rc
  out=$("$@" 2>&1); rc=$?
  if [[ "$rc" == "$want" ]]; then ok "$name"; else bad "$name" "期望退出码 ${want}，实际 $rc
       $(printf '%s' "$out" | tail -3 | tr '\n' ' ')"; fi
  LAST_OUT="$out"
}

# 断言上一次 expect 的输出里包含某个字符串
expect_has() {
  local name="$1" pat="$2"
  if printf '%s' "$LAST_OUT" | grep -q -- "$pat"; then ok "$name"
  else bad "$name" "输出里没有 '$pat'"; fi
}

expect_lacks() {
  local name="$1" pat="$2"
  if printf '%s' "$LAST_OUT" | grep -q -- "$pat"; then bad "$name" "输出里不该出现 '$pat'"
  else ok "$name"; fi
}

R1="$T/sub20M.SRR13849313_1.fastq.gz"
R2="$T/sub20M.SRR13849313_2.fastq.gz"
BL="$T/hg38-blacklist.v2.bed"
OUT="$T/results"

# 给 ⑨ 一个假黑名单，否则它会走"原样拷贝"分支，测不到 intersect
printf 'chr1\t10000\t10100\n' > "$BL"

echo "== 基本接口 =="
expect "--version"            0 "$PIPE" --version
expect_has "--version 有版本号" "atac "
expect "--help"               0 "$PIPE" --help
expect_has "--help 有用法"     "ATAC-seq 双端 FASTQ"
expect "--list-stages"        0 "$PIPE" --list-stages
expect_has "列出了 13 个阶段里的 coverage" "coverage"
expect_has "列出了产物目录"     "06_peaks"
expect "无参数时报错"          1 "$PIPE"
expect "未知参数报错"          2 "$PIPE" --bogus
expect "未知阶段报错"          1 "$PIPE" -1 "$R1" -2 "$R2" --stages foo

echo
echo "== dry-run 全链路 =="
expect "dry-run 退出码 0"      0 "$PIPE" -1 "$R1" -2 "$R2" -o "$OUT" --blacklist "$BL" --dry-run
expect_has "样本名从文件名推导" "样本 sub20M.SRR13849313"
expect_has "① FastQC"          "fastqc -t"
expect_has "② fastp"           "--detect_adapter_for_pe"
expect_has "③ bowtie2 局部比对" "--very-fast-local"
expect_has "④ 过滤表达式用 && 连接" 'rname != "chrM" && rname'
expect_has "④ 排除 chrM"       "chrM"
expect_has "⑤ Tn5 偏移"        "--ATACshift"
expect_has "⑥ 重排"            "srt2.bam"
expect_has "⑦ 分箱阈值 147"     "--nfr-max 147"
expect_has "⑧ MACS2 BAMPE"     "-f BAMPE"
expect_has "⑧ 不传 --shift"     "--nomodel"
expect_lacks "⑧ 没有 --extsize" "--extsize"
expect_has "⑨ 黑名单用 -v"      "intersect -v"
expect_has "⑩ FRiP 双算法"      "--bedtools"
expect_has "⑪ CPM 归一"        "--normalizeUsing CPM"
expect_has "⑫ HOMER 注释"      "annotatePeaks.pl"
expect_has "⑬ motif 用 -mknown" "-mknown"
expect_has "⑭ ATACseqQC"       "atacqc.R"
expect_has "⑮ IGV batch"       "igv_batch.txt"
expect_has "报告阶段"           "make_report.py"

echo
echo "== dry-run 的副作用 =="
if compgen -G "$OUT/.status/*.done" >/dev/null; then
  bad "dry-run 不应落完成标记" "发现 $(ls "$OUT/.status" | tr '\n' ' ')"
else
  ok "dry-run 没有落完成标记"
fi
if [[ -d "$OUT/06_peaks" && -d "$OUT/99_report" ]]; then
  ok "dry-run 建好了目录骨架"
else
  bad "dry-run 没有建目录骨架"
fi
if [[ -s "$OUT/metrics/01_fastqc.tsv" ]]; then
  bad "dry-run 不应写指标文件"
else
  ok "dry-run 没有写指标文件"
fi

echo
echo "== --stages 只跑指定阶段 =="
expect "--stages peaks 退出码 0" 0 "$PIPE" -1 "$R1" -2 "$R2" -o "$OUT" --stages peaks --dry-run
expect_has "包含 peaks 阶段"     "MACS2"
expect_lacks "不含 FastQC"       "fastqc -t"
expect_lacks "不含 TSS"          "computeMatrix"

echo
echo "== 阶段跳过 =="
mkdir -p "$OUT/.status" && : > "$OUT/.status/peaks.done"
expect "已完成阶段被跳过"        0 "$PIPE" -1 "$R1" -2 "$R2" -o "$OUT" --stages peaks --dry-run
expect_has "日志里说了跳过"      "已完成，跳过"
expect "--force 强制重跑"        0 "$PIPE" -1 "$R1" -2 "$R2" -o "$OUT" --stages peaks --dry-run --force
expect_has "--force 后真的跑了"  "MACS2"

echo
echo "== 参数覆盖 =="
expect "-s 覆盖样本名"           0 "$PIPE" -1 "$R1" -2 "$R2" -o "$OUT" -s MySample --dry-run
expect_has "样本名变成 MySample" "MySample"
expect "--nfr-max 覆盖阈值"      0 "$PIPE" -1 "$R1" -2 "$R2" -o "$OUT" --nfr-max 150 --dry-run
expect_has "阈值传到了脚本"       "--nfr-max 150"

echo
echo "== deploy =="
expect "deploy --help"           0 bash "$DEPLOY" --help
expect_has "列出了 genome 目标"  "genome"
expect "deploy doctor"           0 bash "$DEPLOY" doctor
expect_has "doctor 有汇总行"     "必需工具"
expect "deploy 未知目标报错"     1 bash "$DEPLOY" bogus

echo
printf '== 结果：%d 通过 / %d 失败 ==\n' "$PASS" "$FAIL"
if (( FAIL > 0 )); then
  printf '失败项：%s\n' "${FAILED_NAMES[*]}"
  exit 1
fi
