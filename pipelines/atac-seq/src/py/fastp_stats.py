#!/usr/bin/env python3
"""② 从 fastp 的 JSON 报告里抽取可比对的指标。

fastp 的 html 给人看，json 给机器读。这里只取报告里真正会用到的几个数：
修剪前后 reads 数、Q30、重复率、被切掉接头的 reads 数。

输出 5 列 TSV：step  key  value  unit  benchmark
与 atac_pipeline.sh 里的 metric() 完全一致，collect_metrics.py 直接合并。
"""

from __future__ import annotations

import argparse
import json
import os
import sys


def g(d, *path, default=None):
    """安全取嵌套字段。"""
    cur = d
    for k in path:
        if not isinstance(cur, dict) or k not in cur:
            return default
        cur = cur[k]
    return cur


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--json", required=True, help="fastp 的 -j 输出")
    ap.add_argument("--out", required=True, help="要写入的 metrics tsv（覆盖）")
    ap.add_argument("--step", default="02_fastp")
    a = ap.parse_args()

    if not os.path.exists(a.json):
        print(f"[fastp_stats] 找不到 {a.json}", file=sys.stderr)
        return 1

    with open(a.json) as fh:
        d = json.load(fh)

    before = g(d, "summary", "before_filtering", default={}) or {}
    after = g(d, "summary", "after_filtering", default={}) or {}
    filt = g(d, "filtering_result", default={}) or {}
    dup = g(d, "duplication", "rate", default=None)
    adapter = g(d, "adapter_cutting", "adapter_trimmed_reads", default=None)

    rows = []

    def add(key, value, unit="", bench=""):
        if value is None:
            return
        rows.append((a.step, key, value, unit, bench))

    b_reads = before.get("total_reads")
    a_reads = after.get("total_reads")
    if b_reads:
        add("reads_before_trim", int(b_reads), "reads", "")
        add("read_pairs_before_trim", int(b_reads) // 2, "pairs", "")
    if a_reads:
        add("reads_after_trim", int(a_reads), "reads", "")
    if b_reads and a_reads:
        add("reads_retained", round(a_reads / b_reads, 4), "ratio", ">0.95")

    if before.get("q30_rate") is not None:
        add("q30_before_trim", round(float(before["q30_rate"]), 4), "ratio", "")
    if after.get("q30_rate") is not None:
        add("q30_after_trim", round(float(after["q30_rate"]), 4), "ratio", ">0.85")

    if before.get("gc_content") is not None:
        add("gc_content", round(float(before["gc_content"]), 4), "ratio", "~0.41")

    if dup is not None:
        add("duplication_rate_fastp", round(float(dup), 4), "ratio", "<0.30")

    if adapter is not None:
        add("adapter_trimmed_reads", int(adapter), "reads", "")

    if filt.get("low_quality_reads") is not None:
        add("low_quality_reads", int(filt["low_quality_reads"]), "reads", "")
    if filt.get("too_short_reads") is not None:
        add("too_short_reads", int(filt["too_short_reads"]), "reads", "")

    os.makedirs(os.path.dirname(os.path.abspath(a.out)) or ".", exist_ok=True)
    with open(a.out, "w") as fh:
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")

    print(f"[fastp_stats] 写入 {len(rows)} 个指标 → {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
