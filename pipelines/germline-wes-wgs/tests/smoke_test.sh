#!/usr/bin/env bash
# =============================================================================
#  tests/smoke_test.sh —— 冒烟测试（不需要任何生信工具）
#
#  只验证"脚手架"是否正确：参数解析、路径派生、阶段跳过、dry-run 不落标记。
#  真实数据必须上服务器跑，这里连 bwa 都不需要装。
#
#    bash tests/smoke_test.sh
# =============================================================================
set -uo pipefail

PROJECT_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
PIPE="$PROJECT_ROOT/bin/germline-wes-wgs"
DEPLOY="$PROJECT_ROOT/deploy/deploy.sh"

T=$(mktemp -d "${TMPDIR:-/tmp}/wes-smoke.XXXXXX")
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

# 断言输出里包含某个字符串
expect_has() {
  local name="$1" pat="$2"
  if printf '%s' "$LAST_OUT" | grep -q -- "$pat"; then ok "$name"
  else bad "$name" "输出里没有 '$pat'"; fi
}

# ---------------------------------------------------------------- 假数据
mkdir -p "$T/raw" "$T/refs"
printf '>chr17\nACGTACGTACGTACGT\n' > "$T/refs/ref.fa"
printf 'chr17\t16\t0\t16\t17\n'     > "$T/refs/ref.fa.fai"
: > "$T/refs/ref.dict"
: > "$T/refs/ref.fa.bwt"
printf 'chr17\t1\t16\n' > "$T/refs/exome.bed"
gzip -c < /dev/null > "$T/raw/Case1_R1.fastq.gz"
gzip -c < /dev/null > "$T/raw/Case1_R2.fastq.gz"
printf '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr17\tLN:16\n@RG\tID:x\tSM:Case1\n' > "$T/Case1.mdup.bam"

COMMON=(-r "$T/refs/ref.fa" -b "$T/refs/exome.bed" -t 2)
FASTQ=(-1 "$T/raw/Case1_R1.fastq.gz" -2 "$T/raw/Case1_R2.fastq.gz")

# ---------------------------------------------------------------- 0. 源码规范
printf '\n[0] 源码规范\n'
# bash 3.2（macOS 自带）会把 $VAR 后紧跟的全角标点算进变量名，导致 unbound variable。
# 这类 bug 只在老 bash 上炸，很容易漏，所以固化成检查。
HITS=$(LC_ALL=C awk '
  {
    if ($0 ~ /^[ \t]*#/) next
    line = $0
    while (match(line, /\$[A-Za-z_][A-Za-z0-9_]*[^ -~]/)) {
      print FILENAME ":" FNR ": " substr(line, RSTART, RLENGTH)
      line = substr(line, RSTART + RLENGTH)
    }
  }' "$PROJECT_ROOT"/bin/germline-wes-wgs "$PROJECT_ROOT"/src/*.sh \
     "$PROJECT_ROOT"/deploy/*.sh 2>/dev/null)
if [[ -z "$HITS" ]]; then
  ok "没有 \$VAR 后紧跟全角字符（bash 3.2 兼容）"
else
  bad "没有 \$VAR 后紧跟全角字符" "$(printf '%s' "$HITS" | head -3 | tr '\n' ' ')"
fi

# ---------------------------------------------------------------- 1. CLI
printf '\n[1] 命令行接口\n'
expect "帮助退出码 0"            0 "$PIPE" -h
expect_has "帮助里有 caller 说明" "bcftools"
expect "版本输出"                0 "$PIPE" --version
expect "无输入报错"              1 "$PIPE" "${COMMON[@]}"
expect "非法 caller 报错"        1 "$PIPE" "${FASTQ[@]}" "${COMMON[@]}" -c foo
expect "非法阶段名报错"          1 "$PIPE" "${FASTQ[@]}" "${COMMON[@]}" -c bcftools --stages qc,bogus
expect "阶段清单"                0 "$PIPE" --list-stages
expect_has "阶段清单含 annotate" "annotate"

# ---------------------------------------------------------------- 2. dry-run
printf '\n[2] dry-run（不执行、不落标记）\n'
rm -rf "$T/out1"
expect "全流程三线 dry-run" 0 "$PIPE" "${FASTQ[@]}" "${COMMON[@]}" \
       -o "$T/out1" -c all -a none --dry-run
expect_has "打印了 bwa 命令"      "bwa mem"
expect_has "打印了 markdup 顺序"  "samtools fixmate"
expect_has "打印了 GATK HC"       "HaplotypeCaller"
expect_has "打印了 DeepVariant"   "run_deepvariant"
if [[ -z "$(ls -A "$T/out1/.status" 2>/dev/null)" ]]; then
  ok "dry-run 未写完成标记"
else
  bad "dry-run 未写完成标记" "标记目录非空: $(ls "$T/out1/.status")"
fi

printf '\n[3] 已有 BAM + 注释组合\n'
rm -rf "$T/out2"
expect "BAM + GATK + Funcotator" 0 "$PIPE" --bam "$T/Case1.mdup.bam" "${COMMON[@]}" \
       -o "$T/out2" -c gatk -a funcotator --dry-run
expect_has "跳过 qc/trim/align" "bam,call,annotate,report"
expect_has "Funcotator 命令"     "Funcotator"

rm -rf "$T/out3"
expect "三注释器同时开" 0 "$PIPE" --bam "$T/Case1.mdup.bam" "${COMMON[@]}" \
       -o "$T/out3" -c bcftools -a all --stages annotate --dry-run
expect_has "VEP 命令"      "vep --fork"
expect_has "ANNOVAR 命令"  "table_annovar.pl"

# ---------------------------------------------------------------- 4. 断点续跑
printf '\n[4] 断点续跑\n'
mkdir -p "$T/out1/.status"; touch "$T/out1/.status/qc.done"
expect "已完成阶段跳过" 0 "$PIPE" "${FASTQ[@]}" "${COMMON[@]}" \
       -o "$T/out1" -c bcftools --stages qc --dry-run
expect_has "日志提示跳过" "跳过"
expect "--force 强制重跑" 0 "$PIPE" "${FASTQ[@]}" "${COMMON[@]}" \
       -o "$T/out1" -c bcftools --stages qc --force --dry-run
expect_has "--force 后真的打印命令" "fastqc"

# ---------------------------------------------------------------- 5. 预检
printf '\n[5] 预检\n'
# --check 的退出码取决于本机装了什么，这里只验证它能跑完并如实汇报
# "缺工具时怎么报"必须用空桩测：直接跑的话，本机装了什么就报什么，
# 在工具齐全的服务器上这条永远测不出来（会假红）。
mkdir -p "$T/nobin"
LAST_OUT=$(TOOL_SEARCH_DISABLE=1 EXTRA_PATH=/nonexistent PATH="$T/nobin:/usr/bin:/bin" \
           "$PIPE" "${FASTQ[@]}" "${COMMON[@]}" -o "$T/out4" --check 2>&1) || true
ok "--check 能跑完（不因缺工具中断）"
expect_has "预检打印参数" "样本="
expect_has "如实汇报缺失项" "缺少命令"
if [[ -z "$(ls -A "$T/out4/01_qc" 2>/dev/null)" ]]; then ok "--check 未产生数据"; else bad "--check 未产生数据"; fi

# ---------------------------------------------------------------- 6. deploy
printf '\n[6] deploy（dry-run）\n'
expect "deploy 帮助" 0 "$DEPLOY" -h
expect "deploy status" 0 "$DEPLOY" status
# 指向空目录：本机已经放好参考基因组时 deploy 会直接跳过，就测不到下载逻辑了
expect "deploy --dry-run ref" 0 env GENOME_DIR="$T/emptyref" "$DEPLOY" --dry-run ref
expect_has "打印了下载地址" "storage.googleapis.com"
expect "deploy 未知目标报错" 1 "$DEPLOY" bogus

# ---------------------------------------------------------------- 7. 参考基因组探测
printf '\n[7] 参考基因组自动探测（不看文件名、不需要软链接）\n'
mkdir -p "$T/genome"
printf '>chr17\nACGT\n' > "$T/genome/random_name.fa"                  # 无索引
printf '>chr17\nACGT\n' > "$T/genome/Homo_sapiens_assembly38.fasta"    # 有索引
touch "$T/genome/Homo_sapiens_assembly38.fasta.fai" \
      "$T/genome/Homo_sapiens_assembly38.fasta.bwt" \
      "$T/genome/Homo_sapiens_assembly38.dict"

probe() {
  PROJECT_ROOT="$PROJECT_ROOT" bash -c \
    'source "$1/src/00_common.sh"; discover_ref "$2"' _ "$PROJECT_ROOT" "$1" 2>/dev/null
}

got=$(probe "$T/genome")
if [[ "$got" == *"Homo_sapiens_assembly38.fasta" ]]; then
  ok "优先选带索引的 fasta"
else
  bad "优先选带索引的 fasta" "选中了: ${got:-<空>}"
fi

mkdir -p "$T/genome2"; printf '>chr1\nACGT\n' > "$T/genome2/my_ref.fa"
got=$(probe "$T/genome2")
if [[ "$got" == "$T/genome2/my_ref.fa" ]]; then
  ok "只有裸 fasta 也能选中（无索引）"
else
  bad "只有裸 fasta 也能选中" "选中了: ${got:-<空>}"
fi

mkdir -p "$T/genome3"
got=$(probe "$T/genome3")
if [[ -z "$got" ]]; then ok "空目录返回空（不误报）"; else bad "空目录返回空" "选中了: $got"; fi

# 用户显式指定 REF 时必须优先于自动探测
rm -rf "$T/out5"
LAST_OUT=$("$PIPE" "${FASTQ[@]}" -r "$T/refs/ref.fa" -o "$T/out5" -c bcftools --stages qc --dry-run 2>&1) || true
expect_has "显式 -r 优先生效" "参考: $T/refs/ref.fa"

# ---------------------------------------------------------------- 8. bundle 格式兼容
printf '\n[8] GATK bundle：.vcf 和 .vcf.gz 都要认\n'
mkdir -p "$T/bundle"
for n in dbsnp_146.hg38 hapmap_3.3.hg38 \
         Mills_and_1000G_gold_standard.indels.hg38 \
         1000G_phase1.snps.high_confidence.hg38; do
  printf '##fileformat=VCFv4.2\n' > "$T/bundle/$n.vcf"
  : > "$T/bundle/$n.vcf.idx"
done

probe_bundle() {
  BUNDLE="$T/bundle" bash -c 'source "$1/src/00_common.sh"; bundle_file "$2"' _ "$PROJECT_ROOT" "$1"
}

got=$(probe_bundle dbsnp_146.hg38)
if [[ "$got" == "$T/bundle/dbsnp_146.hg38.vcf" ]]; then
  ok "bundle_file 能取到未压缩 .vcf"
else
  bad "bundle_file 能取到未压缩 .vcf" "拿到: ${got:-<空>}"
fi

# 两种格式都在时优先 .gz
printf '##fileformat=VCFv4.2\n' | gzip -c > "$T/bundle/dbsnp_146.hg38.vcf.gz"
got=$(probe_bundle dbsnp_146.hg38)
if [[ "$got" == "$T/bundle/dbsnp_146.hg38.vcf.gz" ]]; then
  ok "两种格式并存时优先 .vcf.gz"
else
  bad "两种格式并存时优先 .vcf.gz" "拿到: ${got:-<空>}"
fi
rm -f "$T/bundle/dbsnp_146.hg38.vcf.gz"

# 端到端：未压缩 bundle 也要走 BQSR 和 Tranche（而不是被跳过）
rm -rf "$T/out7"
LAST_OUT=$("$PIPE" --bam "$T/Case1.mdup.bam" "${COMMON[@]}" -o "$T/out7" \
           -c gatk -a none --bundle "$T/bundle" --stages call --dry-run 2>&1) || true
expect_has "未压缩 bundle 走 BaseRecalibrator" "BaseRecalibrator"
expect_has "未压缩 bundle 走 Tranche 过滤"     "FilterVariantTranches"

# 压缩格式同样要能走通
mkdir -p "$T/bundle_gz"
for n in dbsnp_146.hg38 hapmap_3.3.hg38 \
         Mills_and_1000G_gold_standard.indels.hg38 \
         1000G_phase1.snps.high_confidence.hg38; do
  printf '##fileformat=VCFv4.2\n' | gzip -c > "$T/bundle_gz/$n.vcf.gz"
done
rm -rf "$T/out8"
LAST_OUT=$("$PIPE" --bam "$T/Case1.mdup.bam" "${COMMON[@]}" -o "$T/out8" \
           -c gatk -a none --bundle "$T/bundle_gz" --stages call --dry-run 2>&1) || true
expect_has "压缩 bundle 同样走 BQSR" "BaseRecalibrator"

# 预检要能点清 4/4
LAST_OUT=$("$PIPE" --bam "$T/Case1.mdup.bam" "${COMMON[@]}" -o "$T/out9" \
           -c gatk --bundle "$T/bundle" --check 2>&1) || true
expect_has "预检报告 4/4 已知位点" "4/4"

# ---------------------------------------------------------------- 9. 软件体检
printf '\n[9] 软件体检：先看已有什么，再决定装什么\n'

expect "doctor 能跑" 0 "$DEPLOY" doctor
expect_has "doctor 打印表头" "工具"
# 同上：空桩里测"缺了会怎么提示"，不依赖本机是否装齐
LAST_OUT=$(TOOL_SEARCH_DISABLE=1 EXTRA_PATH=/nonexistent \
           PATH="$T/nobin:/usr/bin:/bin" "$DEPLOY" doctor 2>&1) || true
expect_has "doctor 给出缺失建议" "deploy.sh env"

# 造一套"已装好"的假工具，验证不会重复安装
mkdir -p "$T/bin"
for t in gatk samtools bwa bcftools fastqc fastp bgzip tabix; do
  printf '#!/bin/sh\nexit 0\n' > "$T/bin/$t"; chmod +x "$T/bin/$t"
done

LAST_OUT=$(TOOL_SEARCH_DISABLE=1 PATH="$T/bin:/usr/bin:/bin" "$DEPLOY" doctor 2>&1) || true
expect_has "必需工具齐全时报 8/8" "8/8"
expect_has "齐全时明确说不需要 conda" "不需要 conda"

LAST_OUT=$(TOOL_SEARCH_DISABLE=1 PATH="$T/bin:/usr/bin:/bin" "$DEPLOY" --dry-run env 2>&1) || true
expect_has "env 遇到齐全环境直接跳过" "跳过 conda"

# 缺一个时只补那一个
rm -f "$T/bin/gatk"
LAST_OUT=$(TOOL_SEARCH_DISABLE=1 PATH="$T/bin:/usr/bin:/bin" "$DEPLOY" --dry-run env 2>&1) || true
expect_has "只缺 gatk 时只装 gatk4" "gatk4"
if printf '%s' "$LAST_OUT" | grep -qE '(^| )bwa( |$)'; then
  bad "只补缺失项（不应重复装 bwa）" "输出里出现了 bwa"
else
  ok "只补缺失项（不应重复装 bwa）"
fi

# status 的"软件"行必须和 doctor 同一口径：有工具就算 OK，
# 不再因为不存在名为 wes 的 conda 环境而报缺失（工具在 /bioapps/bin 里很正常）
printf '#!/bin/sh\nexit 0\n' > "$T/bin/gatk"; chmod +x "$T/bin/gatk"
LAST_OUT=$(TOOL_SEARCH_DISABLE=1 PATH="$T/bin:/usr/bin:/bin" "$DEPLOY" status 2>&1) || true
expect_has "status 报软件就绪" "必需工具就绪"
if printf '%s' "$LAST_OUT" | grep -q "conda 环境"; then
  bad "status 不再要求 conda 环境" "仍出现 'conda 环境'"
else
  ok "status 不再要求 conda 环境"
fi

# 缺工具时 status 要如实报缺，并指到 env
rm -f "$T/bin/gatk"
LAST_OUT=$(TOOL_SEARCH_DISABLE=1 PATH="$T/bin:/usr/bin:/bin" "$DEPLOY" status 2>&1) || true
expect_has "status 缺工具时报缺失" "缺失(1)"
expect_has "status 指到 env 补装" "deploy.sh env"

# 软链接形式的参考基因组：du 默认不跟随软链接，会把 3G 的文件显示成 0
mkdir -p "$T/gsrc" "$T/gdir"
dd if=/dev/zero of="$T/gsrc/big.fa" bs=1024 count=2048 2>/dev/null
ln -sf "$T/gsrc/big.fa" "$T/gdir/my_ref.fa"
printf 'chr17\t8\t0\t8\t9\n' > "$T/gdir/my_ref.fa.fai"
printf 'x\n' > "$T/gdir/my_ref.dict"
printf 'x\n' > "$T/gdir/my_ref.fa.bwt"
LAST_OUT=$(GENOME_DIR="$T/gdir" TOOL_SEARCH_DISABLE=1 PATH="$T/bin:/usr/bin:/bin" "$DEPLOY" status 2>&1) || true
if printf '%s' "$LAST_OUT" | grep -q "my_ref.fa (2.0M)"; then
  ok "软链接参考基因组显示真实大小（不是 0）"
else
  bad "软链接参考基因组显示真实大小" "$(printf '%s' "$LAST_OUT" | grep 参考基因组)"
fi

# ---------------------------------------------------------------- 汇总
printf '\n============================================\n'
printf '通过 %d，失败 %d\n' "$PASS" "$FAIL"
if (( FAIL > 0 )); then
  printf '失败项: %s\n' "${FAILED_NAMES[*]}"
  exit 1
fi
printf '全部通过\n'
