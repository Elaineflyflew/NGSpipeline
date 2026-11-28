#!/usr/bin/env python3
"""⑪ 从 GTF 生成 TSS BED（给 deeptools computeMatrix 用）。

只取 protein_coding 基因，每个基因一个 TSS：
  + 链：TSS = start - 1（BED 是 0-based）
  - 链：TSS = end

必须输出 6 列并带上链方向，否则 computeMatrix --referencePoint TSS 会把
负链基因的上下游搞反，TSS 富集图会出现"两侧高、中间低"的假象。
"""

from __future__ import annotations

import argparse
import gzip
import re
import sys

ATTR = re.compile(r'(\S+)\s+"([^"]*)"')


def open_maybe_gz(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--biotype", default="protein_coding")
    ap.add_argument("--max-genes", type=int, default=0, help="0 表示不限制")
    a = ap.parse_args()

    seen = set()
    n = 0
    with open_maybe_gz(a.gtf) as fh, open(a.out, "w") as out:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue

            attrs = dict(ATTR.findall(f[8]))
            biotype = attrs.get("gene_type") or attrs.get("gene_biotype") or ""
            if a.biotype and biotype != a.biotype:
                continue

            chrom, start, end, strand = f[0], int(f[3]), int(f[4]), f[6]
            if strand == "+":
                tss = start - 1
            elif strand == "-":
                tss = end
            else:
                continue

            name = attrs.get("gene_name") or attrs.get("gene_id") or "NA"
            key = (chrom, tss, strand)
            if key in seen:
                continue
            seen.add(key)

            out.write(f"{chrom}\t{tss}\t{tss + 1}\t{name}\t0\t{strand}\n")
            n += 1
            if a.max_genes and n >= a.max_genes:
                break

    print(f"[tss_prep] 写入 {n} 个 TSS → {a.out}")
    if n == 0:
        print("[tss_prep] 警告：一个 TSS 都没解析出来，检查 GTF 的 feature 列和 biotype",
              file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
