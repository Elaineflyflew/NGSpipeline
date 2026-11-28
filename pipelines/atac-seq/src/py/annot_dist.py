#!/usr/bin/env python3
"""⑫ 解析 HOMER annotatePeaks.pl -annStats 的输出，得到 peak 基因组分布。

HOMER 的 annStats 有两代格式，都得认：

  老格式（HOMER ≤ 4.10）
      Annotation          Total       %
      promoter-TSS        20732       27.96
      ...
      Total               74134       100.00

  新格式（HOMER 4.11+，本仓库在课程服务器上遇到的就是这种）
      Annotation  Number of peaks  Total size (bp)  Log2 Ratio  LogP enrichment
      3UTR        840.0            27132416         0.307       -20.538
      ...
      而且文件里是**多段**：先是标准注释（Promoter/Exon/Intron/Intergenic…），
      后面还有重复序列（LINE/SINE/LTR…）、CpG 岛等，每段都以 Annotation 开头。

所以解析逻辑是：以 "Annotation" 开头的行分段，每段内部按名字+计数读；
百分比自己算（新格式根本不写百分比）。第一段是标准注释，写进主 TSV；
其余段写进 <out>.blocks.tsv 备查。

ATAC 的预期：绝大部分 peak 落在 promoter / intron / distal intergenic。
promoter 占比异常高（>60%）要怀疑 peak 被 TSS 附近的高信号带偏。
"""

from __future__ import annotations

import argparse
import os
import sys


def to_num(s: str):
    """'840.0' / '1,234' / '27.96' → float；不是数字返回 None。"""
    s = s.strip().replace(",", "")
    try:
        return float(s)
    except ValueError:
        return None


def parse_blocks(path: str):
    """返回 [(rows, total), ...]；rows 是 [(name, count, size, log2r), ...]。

    total 优先用 HOMER 自己写的 Total 行（老格式有），没有就事后求和（新格式）。
    """
    blocks = []
    cur = []
    cur_total = None
    with open(path, errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            f = line.split("\t")
            if len(f) < 2:
                continue

            head = f[0].strip().lower()
            if head.startswith("annotation"):
                # 新的一段开始
                if cur:
                    blocks.append((cur, cur_total))
                    cur = []
                    cur_total = None
                continue

            count = to_num(f[1])
            if count is None:
                continue
            if head == "total":          # 合计行：留作分母，本身不算一类
                cur_total = count
                continue
            size = to_num(f[2]) if len(f) > 2 else None
            log2r = to_num(f[3]) if len(f) > 3 else None
            cur.append((f[0].strip(), count, size, log2r))

    if cur:
        blocks.append((cur, cur_total))
    return blocks


def plot(rows, out_png, title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = sorted(rows, key=lambda r: r[1], reverse=True)[:12][::-1]
    labels = [r[0] for r in rows]
    pcts = [r[2] for r in rows]

    fig, ax = plt.subplots(figsize=(6.4, max(2.4, 0.32 * len(rows) + 1.2)), dpi=150)
    bars = ax.barh(labels, pcts, color="#3d6b8c")
    for b, p in zip(bars, pcts):
        ax.text(b.get_width() + 0.6, b.get_y() + b.get_height() / 2,
                f"{p:.1f}%", va="center", fontsize=7)
    ax.set_xlabel("% of peaks")
    ax.set_title(title, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--annstats", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--png", default="")
    a = ap.parse_args()

    if not os.path.exists(a.annstats):
        print(f"[annot_dist] 找不到 {a.annstats}", file=sys.stderr)
        return 1

    blocks = parse_blocks(a.annstats)
    if not blocks or not blocks[0][0]:
        print("[annot_dist] 没解析出分布，检查 annStats 文件格式", file=sys.stderr)
        return 1

    # 第一段 = 标准注释（Promoter / Exon / Intron / Intergenic …）
    main_rows, main_total = blocks[0]
    total = main_total if main_total else sum(c for _, c, _, _ in main_rows)
    total = total or 1.0
    rows = [(name, int(c), round(100.0 * c / total, 2)) for name, c, _, _ in main_rows]

    with open(a.out, "w") as fh:
        fh.write("annotation\tcount\tpercent\n")
        for name, count, pct in rows:
            fh.write(f"{name}\t{count}\t{pct}\n")
        # Total 行是给 stage_annot 取 "peaks_annotated" 用的
        fh.write(f"Total\t{int(total)}\t100.00\n")

    # 其余段（重复序列 / CpG 等）单独存，不干扰主表
    if len(blocks) > 1:
        alt = a.out[:-4] + ".blocks.tsv" if a.out.endswith(".tsv") else a.out + ".blocks.tsv"
        with open(alt, "w") as fh:
            fh.write("block\tannotation\tcount\tpercent\ttotal_size\tlog2_ratio\n")
            for bi, (blk, blk_total) in enumerate(blocks):
                t = blk_total if blk_total else sum(c for _, c, _, _ in blk)
                t = t or 1.0
                for name, c, size, l2 in blk:
                    fh.write("%d\t%s\t%d\t%.2f\t%s\t%s\n" % (
                        bi, name, int(c), 100.0 * c / t,
                        "" if size is None else int(size),
                        "" if l2 is None else l2))
        print(f"[annot_dist] 另 {len(blocks) - 1} 段（重复序列等）→ {alt}")

    if a.png:
        try:
            plot(rows, a.png, "Peak distribution (HOMER annStats)")
        except Exception as exc:   # 没装 matplotlib 不该让主链失败
            print(f"[annot_dist] 跳过画图: {exc}", file=sys.stderr)

    top = max(rows, key=lambda r: r[1])
    print(f"[annot_dist] {len(rows)} 类；最大类 {top[0]} {top[2]:.1f}% → {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
