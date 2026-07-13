#!/usr/bin/env python3
# ============================================================
# Fig 5H — plotting (stacked/grouped bars, pixel co-localisation)
# Study      : Distinct CAF lineages and an Am80 stromal-reprogramming
#              strategy in extrahepatic cholangiocarcinoma (Meflin/ISLR).
# Repository : eCCA_AM80  (https://github.com/shirakiyukihiro/eCCA_AM80)
# Input = per_image_summary.csv from fig5h_coexpr.py.
# Paths made relative for public release: put inputs in ./data, outputs in ./output.
# ============================================================

"""
plot_fig5h.py -- Fig 5H figures, styled to match the R house style (00_plot_style.R):
theme_classic look, Okabe-Ito per-patient dots (white-edged, semi-transparent),
light category fills with a grey outline. Descriptive: mean +/- s.e.m., patient = unit.

Produces:
  Fig5H_coexpr.pdf/.png          ** primary ** single 100% horizontal stacked bar
                                 (mean composition) with per-patient dots + s.e.m.
                                 error bars at the internal boundaries.
  Fig5H_coexpr_grouped.pdf/.png  alternative: 3 vertical bars (mean +/- s.e.m. + dots)
  Fig5H_pixel_coloc.pdf/.png     segmentation-free pixel co-localization (corroboration)

Usage:  python plot_fig5h.py per_image_summary.csv [out_dir]
"""
import sys, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

infile = sys.argv[1] if len(sys.argv) > 1 else "per_image_summary.csv"
outdir = sys.argv[2] if len(sys.argv) > 2 else "."
os.makedirs(outdir, exist_ok=True)

OKABE = ["#E69F00", "#56B4E9", "#009E73", "#F0E442",
         "#0072B2", "#D55E00", "#CC79A7", "#999999"]
GREY35, GREY20 = "#595959", "#333333"
plt.rcParams.update({
    "font.size": 9, "font.family": "sans-serif",
    "axes.linewidth": 0.8, "axes.edgecolor": "black",
    "xtick.direction": "out", "ytick.direction": "out",
    "xtick.color": "black", "ytick.color": "black",
    "pdf.fonttype": 42, "ps.fonttype": 42,
})

pi = pd.read_csv(infile, encoding="cp932")
patients = sorted(pi["patient"].unique())
pal = {p: OKABE[i % len(OKABE)] for i, p in enumerate(patients)}

PROP = ["pct_Meflin_only", "pct_double", "pct_Thy1_only"]
PLAB = ["Meflin$^+$ only", "double$^+$", "Thy1$^+$ only"]
PFILL = ["#2CA02C", "#F0C000", "#D62728"]


def summarise(cols):
    per_pat = pi.groupby("patient")[cols].mean().loc[patients]
    return per_pat, per_pat.mean()[cols].values, per_pat.sem()[cols].values


def patient_legend(ax, y=-0.34):
    ax.legend(title="Patient", fontsize=6.5, title_fontsize=7, frameon=False,
              loc="upper center", bbox_to_anchor=(0.5, y),
              ncol=len(patients), columnspacing=0.7, handletextpad=0.15)


def stacked_hbar(fname):
    """Single 100% horizontal stacked bar (mean) + per-patient dots & s.e.m. at boundaries."""
    per_pat, mean, sem = summarise(PROP)
    cum = np.concatenate([[0], np.cumsum(mean)])          # mean boundaries
    per_cum = per_pat[PROP].cumsum(axis=1)                # per-patient cumulative
    fig, ax = plt.subplots(figsize=(6.2, 2.7))
    h = 0.42
    for i in range(len(PROP)):
        ax.barh(0, mean[i], left=cum[i], height=h, color=PFILL[i], alpha=0.55,
                edgecolor="white", linewidth=1.2, zorder=1)
        ax.text(cum[i] + mean[i] / 2, 0, f"{PLAB[i]}\n{mean[i]:.1f} ± {sem[i]:.1f}%",
                ha="center", va="center", fontsize=8, zorder=4)
    yoff = h / 2 + 0.22
    rng = np.random.default_rng(0)
    for b in range(len(PROP) - 1):                        # internal boundaries
        bx = cum[b + 1]; bvals = per_cum.iloc[:, b].values
        ax.plot([bx, bx], [-h / 2, yoff], color="grey", lw=0.6, ls=":", zorder=2)
        ax.errorbar(bx, yoff, xerr=bvals.std(ddof=1) / np.sqrt(len(bvals)),
                    fmt="none", ecolor=GREY20, elinewidth=0.9, capsize=3, zorder=3)
        for j, p in enumerate(patients):
            ax.scatter(bvals[j], yoff + (rng.random() - 0.5) * 0.14, s=28,
                       facecolor=pal[p], edgecolor="white", linewidth=0.5,
                       alpha=0.95, zorder=5, label=(p if b == 0 else None))
    ax.set_xlim(0, 100); ax.set_ylim(-0.5, 0.78); ax.set_yticks([])
    ax.set_xlabel("% of marker-positive cells", fontsize=9)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.set_title(f"Meflin/Thy1 co-expression, human dCCA (n = {len(patients)})",
                 fontsize=9, pad=6)
    patient_legend(ax, y=-0.42)
    fig.savefig(os.path.join(outdir, fname + ".pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(outdir, fname + ".png"), dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"  {fname}: " + " | ".join(f"{l} {m:.1f}±{s:.1f}%"
          for l, m, s in zip(PLAB, mean, sem)))


def grouped_bar(cols, labels, fills, ylab, ylim, title, fname, jitter=0.16):
    per_pat, mean, sem = summarise(cols)
    x = np.arange(len(cols))
    fig, ax = plt.subplots(figsize=(3.4, 3.5))
    ax.bar(x, mean, width=0.5, color=fills, alpha=0.35,
           edgecolor=GREY35, linewidth=0.4, zorder=1)
    ax.errorbar(x, mean, yerr=sem, fmt="none", ecolor=GREY20,
                elinewidth=0.7, capsize=2.5, capthick=0.7, zorder=2)
    rng = np.random.default_rng(0)
    for p in patients:
        vals = per_pat.loc[p, cols].values
        jit = (rng.random(len(cols)) - 0.5) * 2 * jitter
        ax.scatter(x + jit, vals, s=30, facecolor=pal[p], edgecolor="white",
                   linewidth=0.5, alpha=0.85, zorder=3, label=p)
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel(ylab, fontsize=9); ax.set_ylim(*ylim)
    ax.set_title(title, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    patient_legend(ax, y=-0.22)
    fig.savefig(os.path.join(outdir, fname + ".pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(outdir, fname + ".png"), dpi=220, bbox_inches="tight")
    plt.close(fig)


# ---- primary: single 100% horizontal stacked bar ----
stacked_hbar("Fig5H_coexpr")
# ---- alternative: grouped vertical bars ----
grouped_bar(PROP, ["Meflin$^+$\nonly", "double$^+$", "Thy1$^+$\nonly"], PFILL,
            "% of marker-positive cells", (0, 60),
            f"Meflin/Thy1 co-expression\nhuman dCCA (n = {len(patients)})",
            "Fig5H_coexpr_grouped")
# ---- corroboration: pixel co-localization ----
clc = [c for c in ["coloc_jaccard", "manders_M1_MefInThy", "manders_M2_ThyInMef"]
       if c in pi.columns]
if clc:
    grouped_bar(clc, ["Jaccard\noverlap", "Manders\nM1", "Manders\nM2"][:len(clc)],
                ["#7570B3", "#1B9E77", "#D95F02"][:len(clc)],
                "coefficient", (0, 0.6),
                f"Pixel co-localization\n(threshold-free, n = {len(patients)})",
                "Fig5H_pixel_coloc")
print("done.")
