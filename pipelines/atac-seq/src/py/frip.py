#!/usr/bin/env python3
"""⑩ FRiP（Fraction of Reads in Peaks）双算法交叉验证。

同一个数用两条互不依赖的路径算：
  A. samtools view -c -L peaks.bed       —— htslib 的区间索引，快
  B. bedtools bamtobed | intersect -u    —— 把比对转成区间再求交，慢但实现完全不同

两个数应该几乎一致。差得多就说明有一边出了问题（比如 BED 排序、染色体命名、
peak 文件和 BAM 不是同一版本），这时候不该继续往下解读 FRiP。

经验阈值（ENCODE ATAC-seq 标准）：
  > 0.30  优秀
  0.20-0.30  合格
  0.10-0.20  偏低，通常是 signal-to-noise 不够
  < 0.10  基本失败
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys


def sh(args: list[str]) -> str:
    return subprocess.check_output(args, text=True, stderr=subprocess.DEVNULL).strip()


def count_bedtools(bam: str, peaks: str) -> int:
    """bedtools bamtobed → intersect -u → 数行数（用 Python 数，避免 wc 的管道开销）。"""
    p1 = subprocess.Popen(
        ["bedtools", "bamtobed", "-i", bam],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True,
    )
    p2 = subprocess.Popen(
        ["bedtools", "intersect", "-u", "-a", "stdin", "-b", peaks],
        stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True,
    )
    assert p1.stdout is not None
    p1.stdout.close()
    n = 0
    assert p2.stdout is not None
    for _ in p2.stdout:
        n += 1
    p2.wait()
    p1.wait()
    return n


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--peaks", required=True)
    ap.add_argument("--label", required=True, help="all / nfr / mono")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--bedtools", action="store_true", help="额外跑一遍 bedtools 算法")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--step", default="", help="给定时输出 <step>\\t<label>_<key>\\t... 的标准指标格式")
    a = ap.parse_args()

    os.makedirs(a.outdir, exist_ok=True)

    total = int(sh(["samtools", "view", "-c", "-F", "4", a.bam]))
    in_peak = int(sh(["samtools", "view", "-c", "-F", "4", "-L", a.peaks, a.bam]))

    rows = [("total_reads", total, "reads", "")]
    rows.append(("reads_in_peaks", in_peak, "reads", ""))
    frip_a = (in_peak / total) if total else 0.0
    rows.append(("frip_samtools", round(frip_a, 4), "ratio", ">0.20"))

    if a.bedtools:
        try:
            n_b = count_bedtools(a.bam, a.peaks)
            frip_b = (n_b / total) if total else 0.0
            rows.append(("reads_in_peaks_bedtools", n_b, "reads", ""))
            rows.append(("frip_bedtools", round(frip_b, 4), "ratio", ">0.20"))
            rows.append(("frip_delta", round(abs(frip_a - frip_b), 4), "ratio", "<0.02"))
        except FileNotFoundError:
            print("[frip] 没有 bedtools，跳过算法 B", file=sys.stderr)
        except subprocess.CalledProcessError as exc:
            print(f"[frip] bedtools 失败: {exc}", file=sys.stderr)

    out = os.path.join(a.outdir, f"{a.label}.frip.tsv")
    with open(out, "w") as fh:
        for key, val, unit, bench in rows:
            if a.step:
                fh.write(f"{a.step}\t{a.label}_{key}\t{val}\t{unit}\t{bench}\n")
            else:
                fh.write(f"{a.label}\t{key}\t{val}\t{unit}\t{bench}\n")

    print(f"[frip] {a.label}: {in_peak}/{total} = {frip_a:.4f} → {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
