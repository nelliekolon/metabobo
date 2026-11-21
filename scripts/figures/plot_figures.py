#!/usr/bin/env python3
"""Draw the write-up figures from the tracked figure data.

usage: python3 plot_figures.py <repo root> <out dir>
  batch_drift.png  held-out pooled QCs in MTBLS2469, by run day, before and after batch correction
  rt_map.png       retention-time maps from each study onto the reference (anchors used in the fit)
"""
import os, sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

INK, MUTED, BEFORE, AFTER = "#222222", "#777777", "#b0b0b0", "#2f6db5"
STUDIES = {"MSV000084556": "#2f6db5", "MSV000094097": "#d1803a", "MTBLS2469": "#4a9a5b"}


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(colors=MUTED, labelsize=9)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)


def batch_drift(root, out):
    d = pd.read_csv(os.path.join(root, "processed/figures/batch_drift_MTBLS2469.tsv"), sep="\t")
    days = sorted(d.batch.unique())
    fig, ax = plt.subplots(figsize=(6.4, 3.4), dpi=200)
    for stage, color, dx in (("before", BEFORE, -0.12), ("after", AFTER, 0.12)):
        s = d[d.stage == stage]
        x = [days.index(b) + dx for b in s.batch]
        ax.scatter(x, s.deviation_log2, s=22, color=color, label=f"{stage} correction", zorder=3)
    ax.axhline(0, color=MUTED, lw=0.8, zorder=1)
    ax.set_xticks(range(len(days)), [b[5:] for b in days], rotation=0)
    ax.set_xlabel("run day (2020)", color=INK, fontsize=9)
    ax.set_ylabel("pool deviation (log2)", color=INK, fontsize=9)
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    style(ax)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "batch_drift.png"))


def rt_map(root, out):
    fig, ax = plt.subplots(figsize=(4.6, 4.2), dpi=200)
    for study, color in STUDIES.items():
        m = pd.read_csv(os.path.join(root, f"processed/crossstudy/rt_map_{study}.tsv"), sep="\t")
        m = m[m.used.astype(str) == "True"]
        ax.scatter(m.rt_q / 60, m.rt_ref / 60, s=8, color=color, alpha=0.7, label=f"{study} ({len(m)})")
    ax.plot([0, 13], [0, 13], color=MUTED, lw=0.8, ls="--")
    ax.set_xlabel("retention time in study (min)", color=INK, fontsize=9)
    ax.set_ylabel("retention time in reference (min)", color=INK, fontsize=9)
    ax.legend(frameon=False, fontsize=8, title="anchors used", title_fontsize=8)
    style(ax)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "rt_map.png"))


if __name__ == "__main__":
    root, out = sys.argv[1:3]
    os.makedirs(out, exist_ok=True)
    batch_drift(root, out)
    rt_map(root, out)
