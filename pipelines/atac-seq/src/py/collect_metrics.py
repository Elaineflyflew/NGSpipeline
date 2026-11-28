#!/usr/bin/env python3
"""把 metrics/*.tsv 合并成一张 summary.tsv，并按阈值打 PASS / WARN / FAIL。

为什么要有这一步：
    15 个步骤各自写了指标，但"指标"本身不说话。判读需要三样东西放在一起——
    数值、参考阈值、结论。这张表就是中间那层：把散落的数值按同一把尺子量一遍，
    报告里只引用它的 status，不再重复写判断逻辑。
"""

from __future__ import annotations

import argparse
import glob
import os
import re
import sys

# 各步骤的展示顺序（和作业图的编号对齐）
STEP_ORDER = [
    "01_fastqc", "02_fastp", "03_bowtie2", "04_bam_filter", "05_shift",
    "06_resort", "07_fragbin", "08_macs2", "09_blacklist", "10_frip",
    "11_coverage", "12_annotate", "13_motif", "14_atacqc", "15_igv", "99_report",
]

STEP_TITLE = {
    "01_fastqc": "① FastQC 原始质检",
    "02_fastp": "② fastp 修剪",
    "03_bowtie2": "③ bowtie2 比对",
    "04_bam_filter": "④ 过滤 + 排序",
    "05_shift": "⑤ Tn5 校正",
    "06_resort": "⑥ 重排",
    "07_fragbin": "⑦ 片段分箱",
    "08_macs2": "⑧ MACS2",
    "09_blacklist": "⑨ 黑名单",
    "10_frip": "⑩ FRiP",
    "11_coverage": "⑪ coverage + TSS",
    "12_annotate": "⑫ peak 注释",
    "13_motif": "⑬ motif 富集",
    "14_atacqc": "⑭ ATACseqQC",
    "15_igv": "⑮ IGV",
    "99_report": "汇总",
}


def to_number(s: str):
    """把 20k / 1.2e6 / 53.88% 之类转成 float。"""
    s = s.strip()
    m = re.fullmatch(r"(-?\d+(?:\.\d+)?)([kKmM])?", s)
    if m:
        v = float(m.group(1))
        mult = {"k": 1e3, "m": 1e6}.get((m.group(2) or "").lower(), 1)
        return v * mult
    try:
        return float(s)
    except ValueError:
        return None


def judge(value: str, bench: str) -> str:
    """按 benchmark 字符串判 PASS/WARN/FAIL；无法判读就返回 '-'。"""
    if not bench or bench in ("-", "NA"):
        return "-"
    if value in ("NA", "SKIP", ""):
        return "-"
    v = to_number(value)
    if v is None:
        return "-"

    b = bench.strip()
    try:
        if b.startswith(">"):
            return "PASS" if v > float(b[1:]) else "WARN"
        if b.startswith("<"):
            return "PASS" if v < float(b[1:]) else "WARN"
        if b.startswith("~"):
            target = float(b[1:])
            return "PASS" if abs(v - target) <= 0.25 * abs(target) else "WARN"
        if "-" in b and not b.startswith("-"):
            lo_s, hi_s = b.split("-", 1)
            lo, hi = to_number(lo_s), to_number(hi_s)
            if lo is None or hi is None:
                return "-"
            return "PASS" if lo <= v <= hi else "WARN"
    except ValueError:
        return "-"
    return "-"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--metrics", required=True, help="metrics 目录")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    files = sorted(glob.glob(os.path.join(a.metrics, "*.tsv")))
    if not files:
        print(f"[collect_metrics] {a.metrics} 下没有 tsv", file=sys.stderr)
        return 1

    rows = []
    # (step, key) -> 行号：同一个指标重跑过多次时以最后一行为准。
    # 为什么不靠写的时候去重：老版本的 metric() 是纯追加，历史文件里已经
    # 留了重复行（实测 04_bam_filter / 07_fragbin 各重复一遍），这里兜底。
    seen = {}
    for path in files:
        # 文件名可以是 <step>.tsv（单样本）或 <step>.<sample>[.<frac>].tsv（多样本）
        step = os.path.basename(path).split(".")[0]
        with open(path, errors="replace") as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 3:
                    continue
                if f[0] == "step":
                    continue
                _, key, value = f[0], f[1], f[2]
                unit = f[3] if len(f) > 3 else ""
                bench = f[4] if len(f) > 4 else ""
                rec = (step, key, value, unit, bench, judge(value, bench))
                if (step, key) in seen:
                    rows[seen[(step, key)]] = rec
                else:
                    seen[(step, key)] = len(rows)
                    rows.append(rec)

    def sort_key(r):
        idx = STEP_ORDER.index(r[0]) if r[0] in STEP_ORDER else 99
        return (idx, r[1])

    rows.sort(key=sort_key)

    with open(a.out, "w") as fh:
        fh.write("step\tmetric\tvalue\tunit\tbenchmark\tstatus\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")

    n_pass = sum(1 for r in rows if r[5] == "PASS")
    n_warn = sum(1 for r in rows if r[5] == "WARN")
    print(f"[collect_metrics] {len(rows)} 条指标（PASS {n_pass} / WARN {n_warn}）→ {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
