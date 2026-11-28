#!/usr/bin/env python3
"""⑪ 从 computeMatrix 的矩阵里算 TSS enrichment score。

定义（与 ENCODE / ATACseqQC 的 TSSE 口径一致）：
    score = mean(signal in TSS ± window) / mean(signal in 两端 flank 区)

为什么不用 plotProfile 的图去目测：
    图能看出"有没有富集"，但报告里需要一个数才能和阈值比、才能跨样本比。
    这个数和 ATACseqQC 的 TSSEscore 是同一个定义，⑭ 可以拿来做交叉验证。

computeMatrix 输出格式：首行是 @ 开头的 JSON（含 sample_boundaries / bin_size），
后面每行是 chr start end v1 v2 ...。用 sample_boundaries 切列，比硬编码列号稳。
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import sys


def open_maybe_gz(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--matrix", required=True, help="computeMatrix -o 的输出（.gz）")
    ap.add_argument("--flank", type=int, default=2000)
    ap.add_argument("--window", type=int, default=50, help="TSS ±window 作为信号区")
    ap.add_argument("--out", required=True)
    ap.add_argument("--step", default="", help="给定时输出标准 5 列指标")
    a = ap.parse_args()

    bounds = None
    bin_size = None
    n_bins = None
    n_bins_total = None
    region_sum = None
    n_regions = 0

    def scalar(v):
        """deeptools 3.5 起元数据变成了单元素列表（"bin size":[25]），统一取标量。"""
        if isinstance(v, (list, tuple)):
            return v[0] if v else None
        return v

    with open_maybe_gz(a.matrix) as fh:
        for line in fh:
            if line.startswith("@"):
                payload = line[1:].strip()
                if not payload.startswith("{"):
                    continue
                try:
                    meta = json.loads(payload)
                except json.JSONDecodeError:
                    continue
                # sample_boundaries 有 [[0, 160]] 和 [0, 160, 320] 两种写法，都兼容
                sb = meta.get("sample_boundaries")
                if isinstance(sb, list) and sb:
                    if isinstance(sb[0], (list, tuple)):
                        sb = sb[0]
                    try:
                        bounds = [int(x) for x in sb]
                        n_bins_total = bounds[-1]
                    except (TypeError, ValueError):
                        bounds = None
                bs = scalar(meta.get("bin size", meta.get("bin_size")))
                if bs is not None:
                    try:
                        bin_size = int(bs)
                    except (TypeError, ValueError):
                        pass
                continue

            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            # computeMatrix 的行 = 若干注释列（chr/start/end[/name/score/strand]）+ 数据列。
            # 数据列一定在最后，所以先用「总列数 − bin 数」定位；不行再退回
            # sample_boundaries 给的下标（老版本里那是绝对列号）。
            cand = []
            if n_bins_total and len(f) > n_bins_total:
                cand.append(f[len(f) - n_bins_total:])
            if bounds and 0 <= bounds[0] < bounds[-1] <= len(f):
                cand.append(f[bounds[0]:bounds[-1]])
            cand.append(f[6:])

            vals = None
            for raw in cand:
                try:
                    v = [float(x) for x in raw]
                except ValueError:
                    continue
                if len(v) >= 10:
                    vals = v
                    break
            if not vals:
                continue
            if not vals:
                continue
            if region_sum is None:
                region_sum = [0.0] * len(vals)
                n_bins = len(vals)
            if len(vals) != n_bins:
                continue
            for i, v in enumerate(vals):
                region_sum[i] += (0.0 if math.isnan(v) else v)
            n_regions += 1

    if not n_regions or region_sum is None:
        print("[tss_score] 矩阵里没读到数据行", file=sys.stderr)
        return 1

    profile = [s / n_regions for s in region_sum]

    if bin_size is None:
        bin_size = max(1, (2 * a.flank) // n_bins)

    center = n_bins // 2
    half = max(1, int(round(a.window / 2 / bin_size)))
    flank_bins = max(1, int(round(a.flank / 4 / bin_size)))   # 取两端各 1/4 作为背景

    sig = profile[max(0, center - half): min(n_bins, center + half + 1)]
    bg = profile[:flank_bins] + profile[-flank_bins:]

    sig_mean = sum(sig) / len(sig) if sig else 0.0
    bg_mean = sum(bg) / len(bg) if bg else 0.0
    score = (sig_mean / bg_mean) if bg_mean > 0 else float("nan")

    with open(a.out, "w") as out:
        if a.step:
            out.write(f"{a.step}\ttss_enrichment_score\t{score:.3f}\tfold\t>6\n")
            out.write(f"{a.step}\ttss_signal_mean\t{sig_mean:.4f}\t\t\n")
            out.write(f"{a.step}\tflank_signal_mean\t{bg_mean:.4f}\t\t\n")
            out.write(f"{a.step}\tn_tss_regions\t{n_regions}\t\t\n")
        else:
            out.write(f"tss_enrichment_score\t{score:.3f}\n")
            out.write(f"tss_signal_mean\t{sig_mean:.4f}\n")
            out.write(f"flank_signal_mean\t{bg_mean:.4f}\n")
            out.write(f"n_tss_regions\t{n_regions}\n")
            out.write(f"bin_size\t{bin_size}\n")

    print(f"[tss_score] TSS enrichment = {score:.2f}（{n_regions} 个 TSS，bin {bin_size} bp）")
    return 0


if __name__ == "__main__":
    sys.exit(main())
