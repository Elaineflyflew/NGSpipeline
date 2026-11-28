#!/usr/bin/env python3
"""把 summary.tsv + 各步骤产物，拼成一份自包含的 HTML 报告。

这份报告要回答的问题不是"跑完了吗"，而是"跑对了吗"：
  · 每一步的产物在哪（可点开检查）
  · 这一步的关键指标是多少、和阈值比是什么结论
  · 整条链合起来，这个样本能不能用

图直接 base64 嵌进 HTML，单文件就能发给别人看，不依赖外部目录。
"""

from __future__ import annotations

import argparse
import base64
import datetime as dt
import glob
import html
import os
import sys

STEP_TITLE = {
    "01_fastqc": "① FastQC 原始质检",
    "02_fastp": "② fastp 接头修剪 + 质量过滤",
    "03_bowtie2": "③ bowtie2 比对",
    "04_bam_filter": "④ MAPQ 过滤 + 去 chrM + 排序",
    "05_shift": "⑤ alignmentSieve Tn5 校正",
    "06_resort": "⑥ 重排 + 建索引",
    "07_fragbin": "⑦ 按 TLEN 分箱",
    "08_macs2": "⑧ MACS2 peak calling",
    "09_blacklist": "⑨ 去 ENCODE 黑名单",
    "10_frip": "⑩ FRiP 双算法交叉验证",
    "11_coverage": "⑪ bigWig + TSS 富集",
    "12_annotate": "⑫ HOMER peak 注释",
    "13_motif": "⑬ HOMER motif 富集",
    "14_atacqc": "⑭ ATACseqQC 片段分布 / TSSE",
    "15_igv": "⑮ MYC 位点 IGV 截图",
    "99_report": "汇总",
}

METRIC_LABEL = {
    "raw_read_pairs": "原始 reads 对数",
    "reads_before_trim": "修剪前 reads",
    "reads_after_trim": "修剪后 reads",
    "reads_retained": "reads 保留率",
    "q30_before_trim": "修剪前 Q30",
    "q30_after_trim": "修剪后 Q30",
    "duplication_rate_fastp": "fastp 估计重复率",
    "adapter_trimmed_reads": "切除接头的 reads",
    "overall_alignment_rate": "总体比对率",
    "aligned_reads": "比对上 reads",
    "chrm_reads": "chrM reads",
    "chrm_fraction": "chrM 占比",
    "reads_after_filter": "过滤后 reads",
    "duplicate_fraction": "重复率",
    "insert_size_mean": "平均插入长度",
    "record_fields": "BAM 记录字段数",
    "nfr_fraction": "NFR 占比",
    "mono_fraction": "单核小体占比",
    "poly_fraction": "多核小体占比",
    "peaks_all": "peak 数（全部片段）",
    "peaks_nfr": "peak 数（NFR）",
    "peaks_mono": "peak 数（单核小体）",
    "frip_all": "FRiP（全部片段）",
    "frip_nfr": "FRiP（NFR）",
    "frip_mono": "FRiP（单核小体）",
    "frip_plotEnrichment": "FRiP（plotEnrichment）",
    "tss_enrichment_score": "TSS 富集倍数",
    "top_motif_all": "最强 motif",
    "top_motif_nfr": "最强 motif（NFR）",
    "top_motif_mono": "最强 motif（单核小体）",
    "frip_delta": "两算法 FRiP 差值",
    "frip_nfr": "FRiP（NFR）",
    "frip_mono": "FRiP（单核小体）",
    "reads_in_peaks": "peak 内 reads",
    "reads_in_peaks_bedtools": "peak 内 reads（bedtools）",
    "total_reads": "总 reads",
    "frip_bedtools": "FRiP（bedtools）",
    "frip_samtools": "FRiP（samtools）",
    "peaks_all_after_blacklist": "过滤后 peak 数",
    "peaks_nfr_after_blacklist": "过滤后 peak 数（NFR）",
    "peaks_mono_after_blacklist": "过滤后 peak 数（单核小体）",
    "peaks_all_before_blacklist": "过滤前 peak 数",
    "peaks_nfr_before_blacklist": "过滤前 peak 数（NFR）",
    "peaks_mono_before_blacklist": "过滤前 peak 数（单核小体）",
    "blacklist_removed_all": "黑名单剔除比例",
    "blacklist_removed_nfr": "黑名单剔除比例（NFR）",
    "blacklist_removed_mono": "黑名单剔除比例（单核小体）",
    "nfr_reads": "NFR 片段数",
    "mono_reads": "单核小体片段数",
    "poly_reads": "多核小体片段数",
    "reads_after_shift": "shift 后 reads",
    "properly_paired": "正确配对 reads",
    "duplicate_reads": "重复 reads",
    "aligned_reads": "比对上 reads",
    "chrm_reads": "chrM reads",
    "low_quality_reads": "低质量 reads",
    "too_short_reads": "过短 reads",
    "gc_content": "GC 含量",
    "motifs_tested_all": "通过阈值的 motif 数",
    "motifs_tested_nfr": "通过阈值的 motif 数（NFR）",
    "motifs_tested_mono": "通过阈值的 motif 数（单核小体）",
    "tss_signal_mean": "TSS 区平均信号",
    "flank_signal_mean": "两侧平均信号",
    "n_tss_regions": "TSS 个数",
    "bin_size": "bin 大小",
    "median_fragment_length": "片段长度中位数",
    "tss_enrichment_score_r": "TSSE（ATACseqQC）",
    "nfr_fraction_r": "NFR 占比（ATACseqQC）",
    "atacqc": "ATACseqQC",
    "igv_snapshot": "IGV 截图",
    "annotatePeaks": "HOMER 注释",
    "findMotifs": "HOMER motif",
    "trim_skipped": "跳过修剪",
    "atac_shift": "Tn5 校正",
}

CSS = """
:root{
  --bg:#fbfbfa; --fg:#1f2328; --muted:#6b7280; --line:#e5e7eb;
  --pass:#0f7b3f; --warn:#b45309; --fail:#b91c1c; --accent:#1f4e79;
}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--fg);
  font:14px/1.65 -apple-system,BlinkMacSystemFont,"Segoe UI","Noto Sans SC",sans-serif;}
.wrap{max-width:980px;margin:0 auto;padding:40px 24px 80px}
h1{font-size:26px;margin:0 0 6px;letter-spacing:-.01em}
h2{font-size:17px;margin:40px 0 12px;padding-bottom:6px;border-bottom:1px solid var(--line)}
h3{font-size:14px;margin:0 0 8px}
.sub{color:var(--muted);font-size:13px}
.meta{display:flex;flex-wrap:wrap;gap:6px 22px;margin:14px 0 26px;color:var(--muted);font-size:12.5px}
.meta b{color:var(--fg);font-weight:600}
.verdict{background:#fff;border:1px solid var(--line);border-left:4px solid var(--accent);
  padding:16px 18px;border-radius:6px;margin:18px 0 8px}
.verdict p{margin:0 0 8px}
.verdict p:last-child{margin-bottom:0}
table{width:100%;border-collapse:collapse;background:#fff;font-size:13px;
  border:1px solid var(--line);border-radius:6px;overflow:hidden}
th,td{padding:7px 10px;text-align:left;border-bottom:1px solid var(--line);vertical-align:top}
th{background:#f6f7f9;font-weight:600;font-size:12px;color:var(--muted);
  text-transform:uppercase;letter-spacing:.03em}
tr:last-child td{border-bottom:none}
td.num{font-variant-numeric:tabular-nums;white-space:nowrap}
.st{display:inline-block;padding:1px 7px;border-radius:99px;font-size:11px;font-weight:600}
.st-PASS{background:#e7f5ec;color:var(--pass)}
.st-WARN{background:#fdf2e2;color:var(--warn)}
.st-FAIL{background:#fdeaea;color:var(--fail)}
.st--{background:#f1f2f4;color:var(--muted)}
code{background:#f1f2f4;padding:1px 5px;border-radius:3px;font-size:12px;
  font-family:ui-monospace,SFMono-Regular,Menlo,monospace}
.card{background:#fff;border:1px solid var(--line);border-radius:6px;
  padding:14px 16px;margin:14px 0}
.card .path{font-size:12px;color:var(--muted);margin-top:8px;word-break:break-all}
figure{margin:18px 0;background:#fff;border:1px solid var(--line);border-radius:6px;padding:10px}
figure img{width:100%;display:block;border-radius:4px}
figcaption{font-size:12px;color:var(--muted);margin-top:8px}
.note{font-size:12.5px;color:var(--muted)}
ul.tight{margin:6px 0;padding-left:20px}
ul.tight li{margin:3px 0}
footer{margin-top:50px;padding-top:16px;border-top:1px solid var(--line);
  color:var(--muted);font-size:12px}
"""


def read_summary(path: str):
    rows = []
    with open(path, errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(header):
                f += [""] * (len(header) - len(f))
            rows.append(dict(zip(header, f)))
    return rows


def get(rows, step, metric, default="NA"):
    for r in rows:
        if r["step"] == step and r["metric"] == metric:
            return r["value"]
    return default


def fnum(v, nd=2):
    try:
        return f"{float(v):.{nd}f}"
    except (TypeError, ValueError):
        return "—"


def fmt(value: str, unit: str) -> str:
    """ratio 显示成百分比，大整数加千分位，其余原样。"""
    if unit == "ratio":
        try:
            return f"{float(value) * 100:.2f}%"
        except (TypeError, ValueError):
            return str(value)
    try:
        f = float(value)
    except (TypeError, ValueError):
        return str(value)
    if f == int(f) and abs(f) >= 1000:
        return f"{int(f):,}"
    if f == int(f):
        return str(int(f))
    return f"{f:g}"      # 去掉 "11.420" 这种尾随零


def fmt_bench(bench: str, unit: str) -> str:
    """ratio 的阈值也跟着显示成百分比，避免出现"27.34% 对比 <0.05"这种混搭。"""
    if unit == "ratio" and bench and bench[0] in "<>":
        try:
            return f"{bench[0]}{float(bench[1:]) * 100:g}%"
        except ValueError:
            return bench
    return bench


def pct(v):
    try:
        return f"{float(v) * 100:.1f}%"
    except (TypeError, ValueError):
        return "—"


def badge(status: str) -> str:
    s = status if status in ("PASS", "WARN", "FAIL") else "-"
    return f'<span class="st st-{s}">{s}</span>'


def embed(path: str, caption: str, max_mb: float = 4.0) -> str:
    if not path or not os.path.exists(path):
        return ""
    if os.path.getsize(path) > max_mb * 1024 * 1024:
        return f'<figure><figcaption>{html.escape(caption)}（文件过大，见 {html.escape(path)}）</figcaption></figure>'
    with open(path, "rb") as fh:
        b64 = base64.b64encode(fh.read()).decode()
    return (f'<figure><img alt="{html.escape(caption)}" src="data:image/png;base64,{b64}">'
            f'<figcaption>{html.escape(caption)}</figcaption></figure>')


def build_conclusions(rows) -> list[str]:
    """把关键指标翻成人话。每条都是"数值 + 判据 + 结论"。"""
    out = []

    chrm = get(rows, "03_bowtie2", "chrm_fraction")
    if chrm != "NA":
        v = float(chrm)
        if v < 0.05:
            q = "优秀（<5%）"
        elif v < 0.10:
            q = "良好（5–10%）"
        elif v < 0.20:
            q = "可接受（10–20%）"
        else:
            q = "偏差（>20%）"
        out.append(f"线粒体占比 {pct(chrm)}，属于<b>{q}</b>。"
                   "ATAC 建库时线粒体会被一起打碎，这个比例高会稀释核内信号，"
                   "所以第④步必须去掉 chrM 再 call peak。")

    rate = get(rows, "03_bowtie2", "overall_alignment_rate")
    if rate != "NA":
        out.append(f"bowtie2 总体比对率 {fnum(rate, 2)}%。"
                   "<code>--very-fast-local</code> 允许局部比对，"
                   "ATAC 的 reads 常在接头/低复杂度处截断，这个模式下比对率通常 >95%。")

    nfr = get(rows, "07_fragbin", "nfr_fraction")
    mono = get(rows, "07_fragbin", "mono_fraction")
    poly = get(rows, "07_fragbin", "poly_fraction")
    if nfr != "NA":
        v = float(nfr)
        if v >= 0.30:
            q = "正常，说明酶切充分、开放区被有效富集"
        elif v >= 0.20:
            q = "偏低，可能酶切不足或细胞通透性差"
        else:
            q = "明显偏低，提示建库时转座效率不够"
        out.append(f"片段分布：NFR {pct(nfr)} / 单核小体 {pct(mono)} / 多核小体 {pct(poly)}。"
                   f"NFR 占 {pct(nfr)}，{q}。三个峰（<100、~200、~400 bp）是否分明，"
                   "是判断 ATAC 建库成功与否最直接的证据。")

    frip = get(rows, "10_frip", "frip_all")
    if frip != "NA":
        v = float(frip)
        if v > 0.30:
            q = "优秀（ENCODE 标准 >0.30）"
        elif v > 0.20:
            q = "合格（0.20–0.30）"
        elif v > 0.10:
            q = "偏低，signal-to-noise 不足"
        else:
            q = "失败（<0.10）"
        out.append(f"FRiP = {pct(frip)}，<b>{q}</b>。"
                   "这个数说明“测到的 reads 有多大比例真的落在开放区”，"
                   "是 ATAC 数据可用性的第一道硬门槛。")

    delta = get(rows, "10_frip", "frip_delta")
    if delta != "NA":
        out.append(f"samtools 与 bedtools 两条独立算法算出的 FRiP 相差 {fnum(delta, 4)}，"
                   "两者一致说明 peak 文件与 BAM 的染色体命名、坐标体系是对齐的。"
                   "如果差很多，先查 BED 是否排序、chr 前缀是否一致，而不是急着解释生物学。")

    tss = get(rows, "11_coverage", "tss_enrichment_score")
    if tss != "NA":
        v = float(tss)
        q = ("合格（>6，ENCODE 标准）" if v > 6
             else ("边缘（3–6）" if v > 3 else "不合格（<3）"))
        out.append(f"TSS 富集倍数 {fnum(tss, 2)}×，<b>{q}</b>。"
                   "开放染色质在 TSS 两侧形成明显的“V 形”缺口，"
                   "这个倍数低通常意味着信号被背景拉平了。")

    npk = get(rows, "09_blacklist", "peaks_all_after_blacklist")
    if npk != "NA":
        out.append(f"黑名单过滤后保留 {int(float(npk)):,} 个 peak。"
                   "黑名单区域（着丝粒、端粒、rDNA 等高拷贝区）会产生大量假信号，"
                   "不去掉会让“peak 数很多”变成假象。")

    bl = get(rows, "09_blacklist", "blacklist_removed_all")
    if bl != "NA":
        out.append(f"被黑名单剔除的比例 {fnum(bl, 2)}%。"
                   "正常情况下这个数应该在 1% 以内；偏高说明比对大量堆在重复区。")

    motif = get(rows, "13_motif", "top_motif_all")
    if motif not in ("NA", ""):
        out.append(f"最强富集 motif 是 <b>{html.escape(motif)}</b>。"
                   "ATAC 里 CTCF / NFY / SP1-KLF / ETS 家族是开放染色质的标志性 motif，"
                   "如果 top 列表里出现这些，说明 peak 中心确实落在调控元件上。")

    # ---------- 三指标合看 ----------
    # 任何单一 QC 图都只是"体检一项"：TSS 富集、FRiP、片段周期三者
    # 分别对应"信号在不在调控区""信号集中不集中""建库成不成功"，
    # 缺一项都不能下"这批数据能用"的结论。
    trio = []
    if nfr != "NA":
        trio.append(("片段周期", float(nfr) >= 0.30, "NFR %s" % pct(nfr)))
    if frip != "NA":
        trio.append(("FRiP", float(frip) > 0.20, pct(frip)))
    if tss != "NA":
        trio.append(("TSS 富集", float(tss) > 6, "%s×" % fnum(tss, 2)))

    if len(trio) == 3:
        good = sum(1 for _, ok, _ in trio if ok)
        detail = "；".join("%s %s %s" % (name, val, "达标" if ok else "未达标")
                          for name, ok, val in trio)
        if good == 3:
            verdict = "三项全达标，可以进入生物学解读"
        elif good == 2:
            verdict = ("两项达标、一项偏离——先查偏离那一项的成因，"
                       "不要用另外两项的正常去抵消它")
        else:
            verdict = ("达标项不足两项，这批数据不适合直接做下游解读，"
                       "先回建库/比对环节排查")
        out.append("三指标合看：%s —— <b>%s</b>。" % (detail, verdict))

    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--outdir", required=True, help="99_report 目录")
    ap.add_argument("--title", default="ATAC-seq 分析报告")
    a = ap.parse_args()

    rows = read_summary(a.summary)
    results_root = os.path.dirname(os.path.abspath(a.outdir))
    S = a.sample

    n_pass = sum(1 for r in rows if r["status"] == "PASS")
    n_warn = sum(1 for r in rows if r["status"] == "WARN")

    # ---------- 分步骤表 ----------
    by_step: dict[str, list] = {}
    for r in rows:
        by_step.setdefault(r["step"], []).append(r)

    cards = []
    for step in sorted(by_step, key=lambda s: list(STEP_TITLE).index(s) if s in STEP_TITLE else 99):
        title = STEP_TITLE.get(step, step)
        body = [f'<div class="card"><h3>{html.escape(title)}</h3>']
        body.append('<table><tr><th>指标</th><th>数值</th><th>参考</th><th>判读</th></tr>')
        for r in by_step[step]:
            label = METRIC_LABEL.get(r["metric"], r["metric"])
            unit = f' <span class="sub">{html.escape(r["unit"])}</span>' \
                if r["unit"] and r["unit"] != "ratio" else ""
            body.append(
                f'<tr><td>{html.escape(label)}</td>'
                f'<td class="num">{html.escape(fmt(r["value"], r["unit"]))}{unit}</td>'
                f'<td class="num">{html.escape(fmt_bench(r["benchmark"], r["unit"])) or "—"}</td>'
                f'<td>{badge(r["status"])}</td></tr>'
            )
        body.append("</table>")
        body.append(f'<div class="path">产物目录：<code>results/{step}/</code></div></div>')
        cards.append("".join(body))

    # ---------- 图 ----------
    def find(*patterns):
        for p in patterns:
            hit = sorted(glob.glob(os.path.join(results_root, p)))
            if hit:
                return hit[0]
        return ""

    figs = [
        embed(find(f"05_fragbin/{S}.fragsize_dist.png"),
              "⑦ 片段长度分布：<147 bp（NFR）/ 147–294（单核小体）/ >294（多核小体）"),
        embed(find(f"07_frip/{S}.plotFingerprint.png"),
              "⑩ plotFingerprint：曲线越偏离对角线，signal-to-noise 越好"),
        embed(find(f"07_frip/{S}.plotEnrichment.png"),
              "⑩ plotEnrichment：peak 区域 reads 的累积比例（第二种 FRiP 算法）"),
        embed(find("10_annot/all/annotation_distribution.png"),
              "⑫ peak 基因组分布（HOMER annotatePeaks）"),
        embed(find(f"13_igv/{S}.MYC.png"),
              "⑮ MYC 位点：bigWig 信号 + MACS2 peak 轨道"),
    ]
    figs = [f for f in figs if f]

    conclusions = build_conclusions(rows)
    concl_html = "".join(f"<p>{c}</p>" for c in conclusions) or "<p>指标不足，无法给出结论。</p>"

    now = dt.datetime.now().strftime("%Y-%m-%d %H:%M")
    doc = f"""<!DOCTYPE html>
<html lang="zh-CN"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{html.escape(a.title)} — {html.escape(S)}</title>
<style>{CSS}</style></head><body><div class="wrap">

<h1>{html.escape(a.title)}</h1>
<div class="sub">样本 {html.escape(S)} · SRA SRR13849313（Huh7.5_HCV_ATAC）· hg38 / GENCODE v43</div>
<div class="meta">
  <span>生成时间 <b>{now}</b></span>
  <span>指标总数 <b>{len(rows)}</b></span>
  <span>达标 <b style="color:var(--pass)">{n_pass}</b></span>
  <span>需关注 <b style="color:var(--warn)">{n_warn}</b></span>
</div>

<div class="verdict">
{concl_html}
</div>

<h2>关键判读</h2>
<table>
<tr><th>环节</th><th>指标</th><th>数值</th><th>阈值</th><th>判读</th></tr>
""" + "".join(
        f'<tr><td>{html.escape(STEP_TITLE.get(r["step"], r["step"]))}</td>'
        f'<td>{html.escape(METRIC_LABEL.get(r["metric"], r["metric"]))}</td>'
        f'<td class="num">{html.escape(fmt(r["value"], r["unit"]))}</td>'
        f'<td class="num">{html.escape(fmt_bench(r["benchmark"], r["unit"])) or "—"}</td>'
        f'<td>{badge(r["status"])}</td></tr>'
        for r in rows if r["benchmark"] and r["status"] in ("PASS", "WARN", "FAIL")
    ) + f"""
</table>

<h2>逐步骤产物与指标</h2>
<p class="note">每一步都落一个编号目录：产物证明"跑过了"，指标证明"跑对了"。</p>
{''.join(cards)}

<h2>图</h2>
{''.join(figs) or '<p class="note">没有找到图（可能相关步骤被跳过）。</p>'}

<h2>方法与判据</h2>
<ul class="tight">
  <li><b>比对</b>：bowtie2 <code>--local --very-fast-local</code>，<code>-X 2000</code> 容纳跨核小体片段。</li>
  <li><b>过滤</b>：<code>MAPQ ≥ 20</code> 保留唯一比对，去 chrM 与组装碎片。</li>
  <li><b>Tn5 校正</b>：<code>alignmentSieve --ATACshift</code>（+4/−5），之后必须重排——
      该工具只保留前 11 个字段并打乱坐标顺序。</li>
  <li><b>分箱</b>：按 |TLEN| 分 NFR / 单核小体 / 多核小体；必须用已偏移的 BAM，否则 147 bp 阈值错位 9 bp。</li>
  <li><b>peak</b>：MACS2 <code>-f BAMPE --nomodel</code>，BAMPE 模式下 shift/extsize 会被静默忽略，所以不传。</li>
  <li><b>FRiP</b>：samtools 与 bedtools 两条独立路径，差值应接近 0。</li>
  <li><b>TSS 富集</b>：TSS ±50 bp 均值 ÷ 两端 flank 均值，与 ATACseqQC 的 TSSE 同定义。</li>
</ul>

<footer>
由 <code>atac_pipeline.sh</code> 自动生成 · 指标表 <code>summary.tsv</code> ·
原始 QC 汇总见 <code>multiqc_report.html</code>
</footer>
</div></body></html>
"""

    out = os.path.join(a.outdir, "report.html")
    with open(out, "w") as fh:
        fh.write(doc)
    print(f"[make_report] 报告 → {out}（{len(conclusions)} 条结论，{len(figs)} 张图）")
    return 0


if __name__ == "__main__":
    sys.exit(main())
