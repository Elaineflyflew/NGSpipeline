#!/usr/bin/env python3
"""⑦ 按 TLEN（模板长度）把 BAM 拆成 NFR / 单核小体 / 多核小体三份。

为什么按 TLEN 分：
    ATAC 的双端片段长度直接反映跨越了几个核小体——
      < 147 bp     没被核小体缠住（NFR，nucleosome free region）
      147-294 bp   跨 1 个核小体
      > 294 bp     跨 2 个及以上
    NFR 那一档最富集真正的开放染色质，单核小体档的 peak 反而更像"核小体占位"。

实现上的两个坑（原脚本已踩过，这里固化下来）：
  1. TLEN 有正负号，取决于哪条 mate 在前 → 必须取绝对值
  2. 必须在【已 Tn5 偏移】的 srt2.bam 上分档，否则 147 bp 阈值错位 9 bp

为什么不用 samtools view -e 'abs(tlen)<147'：
   依赖 samtools 的过滤表达式函数，老版本（<1.16）不支持，会静默出错或直接报错。
   这里只用 samtools view -h / -b 这两个任何版本都有的行为，可移植性优先。
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile


def samtools(args: list[str]) -> int:
    return subprocess.call(["samtools"] + args)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bam", required=True, help="已 Tn5 偏移且排好序的 BAM")
    ap.add_argument("--outdir", required=True, help="统计表/图的输出目录")
    ap.add_argument("--bamdir", default="", help="三份 BAM 的输出目录（默认同 --outdir）")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--nfr-max", type=int, default=147)
    ap.add_argument("--mono-max", type=int, default=294)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--hist-max", type=int, default=1000, help="片段分布图画到多长")
    ap.add_argument("--no-plot", action="store_true")
    a = ap.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    bamdir = a.bamdir or a.outdir
    os.makedirs(bamdir, exist_ok=True)

    bins = {
        "nfr":  f"{a.sample}_nucleosome_free.bam",
        "mono": f"{a.sample}_mono_nucleosome.bam",
        "poly": f"{a.sample}_multi_nucleosome.bam",
    }
    paths = {k: os.path.join(bamdir, v) for k, v in bins.items()}

    # 一次遍历：同时写三份 SAM 和长度直方图
    tmpdir = tempfile.mkdtemp(prefix=".fragbin_", dir=a.outdir)
    sams = {k: os.path.join(tmpdir, f"{k}.sam") for k in bins}
    handles = {k: open(p, "w") for k, p in sams.items()}

    hist = [0] * (a.hist_max + 1)
    counts = {"nfr": 0, "mono": 0, "poly": 0, "unpaired": 0}
    tlen_sum = {"nfr": 0, "mono": 0, "poly": 0}

    header_lines = []

    proc = subprocess.Popen(
        ["samtools", "view", "-h", a.bam],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, bufsize=1 << 20,
    )
    assert proc.stdout is not None

    for line in proc.stdout:
        if line.startswith("@"):
            header_lines.append(line)
            continue

        f = line.split("\t", 10)
        if len(f) < 9:
            continue
        flag = int(f[1])
        # 只处理正确配对的 read1，避免同一片段被数两次
        if not (flag & 0x2):
            counts["unpaired"] += 1
            continue
        if not (flag & 0x40):   # 只留 read1：Python 里 & 的优先级低于 ==，必须加括号
            continue

        try:
            tlen = abs(int(f[8]))
        except ValueError:
            continue
        if tlen == 0:
            continue

        if tlen <= a.hist_max:
            hist[tlen] += 1

        if tlen < a.nfr_max:
            key = "nfr"
        elif tlen <= a.mono_max:
            key = "mono"
        else:
            key = "poly"

        counts[key] += 1
        tlen_sum[key] += tlen
        handles[key].write(line)

    proc.wait()
    if proc.returncode != 0:
        print("[frag_bin] samtools view 失败", file=sys.stderr)
        return 1

    # 每个分档：header + records → BAM → 索引
    for key, h in handles.items():
        h.close()
        with open(sams[key]) as body, open(sams[key] + ".full", "w") as out:
            out.writelines(header_lines)
            for line in body:
                out.write(line)
        os.replace(sams[key] + ".full", sams[key])
        if samtools(["view", "-@", str(a.threads), "-b", "-o", paths[key], sams[key]]) != 0:
            print(f"[frag_bin] 生成 {paths[key]} 失败", file=sys.stderr)
            return 1
        if samtools(["index", "-@", str(a.threads), paths[key]]) != 0:
            print(f"[frag_bin] 索引 {paths[key]} 失败", file=sys.stderr)
            return 1
        os.remove(sams[key])
    os.rmdir(tmpdir)

    # 片段长度直方图
    hist_path = os.path.join(a.outdir, f"{a.sample}.fragsize_hist.tsv")
    with open(hist_path, "w") as fh:
        fh.write("fragment_length\tcount\n")
        for i, c in enumerate(hist):
            if c:
                fh.write(f"{i}\t{c}\n")

    # 分档统计
    total = counts["nfr"] + counts["mono"] + counts["poly"]
    stat_path = os.path.join(a.outdir, f"{a.sample}.fragbin_stats.tsv")
    with open(stat_path, "w") as fh:
        fh.write("bin\tfragments\tfraction\tmean_length\tmedian_length\n")
        for key in ("nfr", "mono", "poly"):
            n = counts[key]
            frac = (n / total) if total else 0.0
            mean = (tlen_sum[key] / n) if n else 0.0
            med = median_from_hist(hist, counts, key, a.nfr_max, a.mono_max)
            fh.write(f"{key}\t{n}\t{frac:.4f}\t{mean:.1f}\t{med}\n")

    print(f"[frag_bin] NFR={counts['nfr']} mono={counts['mono']} poly={counts['poly']} "
          f"(未配对跳过 {counts['unpaired']})")

    if not a.no_plot:
        try:
            plot(hist, paths, hist_path.replace("_hist.tsv", "_dist.png"))
        except Exception as exc:  # matplotlib 缺失不该让主链失败
            print(f"[frag_bin] 跳过画图: {exc}", file=sys.stderr)

    return 0


def median_from_hist(hist, counts, key, nfr_max, mono_max):
    """从直方图里求分档中位数：省内存，且对 2500 万条记录足够精确。"""
    lo, hi = (0, nfr_max - 1) if key == "nfr" else \
             (nfr_max, mono_max) if key == "mono" else (mono_max + 1, len(hist) - 1)
    target = counts[key] // 2
    seen = 0
    for i in range(max(lo, 0), min(hi, len(hist) - 1) + 1):
        seen += hist[i]
        if seen >= target:
            return i
    return 0


def plot(hist: list[int], paths: dict, out_png: str) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    xs = list(range(len(hist)))
    ys = hist
    fig, ax = plt.subplots(figsize=(7, 3.2), dpi=150)
    ax.plot(xs, ys, lw=0.8, color="#2c3e50")
    ax.fill_between(xs, ys, color="#5d8aa8", alpha=0.35)

    # 标出两条分界线和 NFR / mono / poly 区间
    for x, label in ((147, "147"), (294, "294")):
        ax.axvline(x, color="#c0392b", ls="--", lw=0.9)
        ax.text(x, max(ys) * 0.95, f" {label} bp", color="#c0392b", fontsize=7, va="top")
    ax.set_xlim(0, min(len(hist) - 1, 800))
    ax.set_xlabel("fragment length (bp)")
    ax.set_ylabel("fragments")
    ax.set_title("Fragment size distribution (NFR / mono- / poly-nucleosome)")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)
    print(f"[frag_bin] 片段分布图 → {out_png}")


if __name__ == "__main__":
    sys.exit(main())
