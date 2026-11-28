#!/usr/bin/env python3
"""⑬ 把 HOMER knownResults.txt 整理成一张能直接看的表。

HOMER 的输出列数随版本变化（有时 9 列、有时多一列 "Motif File"），
所以这里按表头名定位，不硬编码列号。

fold = 目标序列中的占比 / 背景序列中的占比，比单看 p 值更能说明"富集有多强"：
p 值会被序列数放大，fold 不会。ATAC 里常见的高 fold motif 是
CTCF、NFY、SP1/KLF、ETS、AP-1 —— 前三个是开放染色质的标志。
"""

from __future__ import annotations

import argparse
import os
import sys


def find_col(header: list[str], *keys) -> int:
    for i, h in enumerate(header):
        hl = h.strip().lower()
        if all(k in hl for k in keys):
            return i
    return -1


def pct(s: str) -> float:
    s = s.strip().rstrip("%")
    try:
        return float(s)
    except ValueError:
        return float("nan")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--known", required=True, help="HOMER 的 knownResults.txt")
    ap.add_argument("--out", required=True)
    ap.add_argument("--top", type=int, default=15)
    a = ap.parse_args()

    if not os.path.exists(a.known):
        print(f"[motif_top] 找不到 {a.known}", file=sys.stderr)
        return 1

    with open(a.known, errors="replace") as fh:
        lines = [ln.rstrip("\n") for ln in fh if ln.strip()]
    if len(lines) < 2:
        print("[motif_top] knownResults.txt 只有表头，说明没有 motif 通过阈值",
              file=sys.stderr)
        return 1

    header = lines[0].split("\t")
    i_name = find_col(header, "motif", "name")
    i_cons = find_col(header, "consensus")
    i_p = find_col(header, "p-value")
    i_q = find_col(header, "q-value")
    i_tgt = find_col(header, "target", "%")
    i_bg = find_col(header, "background", "%")

    if i_name < 0 or i_tgt < 0 or i_bg < 0:
        print(f"[motif_top] 表头不认识: {header}", file=sys.stderr)
        return 1

    rows = []
    for ln in lines[1:]:
        f = ln.split("\t")
        if len(f) <= max(i_name, i_tgt, i_bg):
            continue
        tgt, bg = pct(f[i_tgt]), pct(f[i_bg])
        fold = (tgt / bg) if bg and bg == bg and bg > 0 else float("nan")
        rows.append({
            "motif": f[i_name].strip(),
            "consensus": f[i_cons].strip() if i_cons >= 0 else "",
            "p": f[i_p].strip() if i_p >= 0 else "",
            "q": f[i_q].strip() if i_q >= 0 else "",
            "target_pct": tgt,
            "background_pct": bg,
            "fold": fold,
        })

    rows.sort(key=lambda r: (-(r["fold"] if r["fold"] == r["fold"] else 0)))

    with open(a.out, "w") as fh:
        fh.write("rank\tmotif\tconsensus\tp_value\tq_value\ttarget_pct\tbackground_pct\tfold\n")
        for i, r in enumerate(rows[:a.top], 1):
            fh.write("\t".join([
                str(i), r["motif"], r["consensus"], r["p"], r["q"],
                f"{r['target_pct']:.2f}", f"{r['background_pct']:.2f}",
                f"{r['fold']:.2f}",
            ]) + "\n")

    print(f"[motif_top] {len(rows)} 个 motif 通过阈值，最高 fold："
          f"{rows[0]['motif']} ({rows[0]['fold']:.1f}×)" if rows else "[motif_top] 无结果")
    return 0


if __name__ == "__main__":
    sys.exit(main())
