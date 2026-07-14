#!/usr/bin/env python
# ============================================================
# Figure 4G — spatial map of tdTomato+ Meflin-lineage cells in the orthotopic tumour
# Study      : Distinct CAF lineages and an Am80 stromal-reprogramming strategy
#              in extrahepatic cholangiocarcinoma (Meflin/ISLR).
# Repository : eCCA_AM80  (https://github.com/shirakiyukihiro/eCCA_AM80)
#
# Renders the published Fig. 4G panel EXACTLY as it appears in the manuscript:
#   - every detected cell centroid as a point (XY projection of the 3D coordinates)
#   - two kernel-density contour lines: 92nd percentile (pink) and 99th percentile (red)
#     of the density values on the evaluation grid, i.e. the contours enclose the
#     densest 8% and the densest 1% of the mapped area, respectively.
#
# NOTE (replaces figure_kde_panelB_style.py):
#   The previous script additionally drew (a) numbered markers at the centroids of the
#   ten largest DBSCAN clusters and (b) an inset schematic of the cluster / focal-cluster
#   / noise classification. Neither appears in the published figure, and the figure legend
#   does not describe them. They are removed here so that the code reproduces the
#   published panel. The DBSCAN / Monte-Carlo cluster statistics are still computed by
#   analyze_tumor_clones_DoG_Zcrop.py and are reported as numbers in the Results text,
#   not as graphical elements of Fig. 4G.
#
# Input      : output/BDC_tumor_DoG_Zcrop_clones.npz  (from analyze_tumor_clones_DoG_Zcrop.py)
# Output     : output/Figure_4G_KDE_map.svg / .png
# ============================================================

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from scipy.stats import gaussian_kde

# ============================================================
# Configuration  (values as reported in the Supplementary Methods)
# ============================================================
NPZ_PATH = 'output/BDC_tumor_DoG_Zcrop_clones.npz'
OUT_DIR = 'output'

KDE_BW = 0.03            # scipy gaussian_kde scalar bandwidth factor
GRID_N = 250             # evaluation grid: GRID_N x GRID_N
GRID_MARGIN_UM = 50      # grid extends this far beyond the data range
CONTOUR_PCTS = [92, 99]  # -> densest 8% and densest 1% of the mapped area
CONTOUR_COLORS = ['#FFB6C1', '#DC143C']  # light pink, crimson
CONTOUR_WIDTHS = [1.2, 2.2]

os.makedirs(OUT_DIR, exist_ok=True)

# ============================================================
# Load detected cell centroids
# ============================================================
d = np.load(NPZ_PATH)
cents = d['centroids_um'].astype(np.float64)   # (N, 3) as (Z, Y, X) in micrometres
n_cells = len(cents)
print(f'Loaded {n_cells:,} cell centroids')

# XY projection: column 2 = X, column 1 = Y
xy = cents[:, [2, 1]]

# ============================================================
# Kernel-density estimate on the XY projection
# ============================================================
kde = gaussian_kde(xy.T, bw_method=KDE_BW)

x_range = (xy[:, 0].min() - GRID_MARGIN_UM, xy[:, 0].max() + GRID_MARGIN_UM)
y_range = (xy[:, 1].min() - GRID_MARGIN_UM, xy[:, 1].max() + GRID_MARGIN_UM)
X, Y = np.meshgrid(np.linspace(*x_range, GRID_N),
                   np.linspace(*y_range, GRID_N))
Z = kde(np.vstack([X.ravel(), Y.ravel()])).reshape(X.shape)

# Iso-density levels: the p-th percentile of the grid density values is the level whose
# contour encloses the densest (100 - p)% of the mapped area.
levels = [np.percentile(Z, p) for p in CONTOUR_PCTS]
for p, lv in zip(CONTOUR_PCTS, levels):
    print(f'  contour: {p}th percentile of grid density '
          f'(densest {100 - p}% of mapped area), level = {lv:.3e}')

# ============================================================
# Figure
# ============================================================
fig, ax = plt.subplots(figsize=(8, 8))

# all detected cells
ax.scatter(xy[:, 0], xy[:, 1], s=3, color='gray', alpha=0.35, zorder=1)

# two-level density contour lines (lines only, not a filled heatmap)
for lv, cc, lw in zip(levels, CONTOUR_COLORS, CONTOUR_WIDTHS):
    ax.contour(X, Y, Z, levels=[lv], colors=[cc], linewidths=lw, zorder=2)

ax.set_xlabel('X (\u00b5m)', fontsize=12)
ax.set_ylabel('Y (\u00b5m)', fontsize=12)
ax.set_aspect('equal')
ax.invert_yaxis()
ax.grid(alpha=0.3)

handles = [
    mlines.Line2D([], [], color='gray', marker='o', markersize=6,
                  linestyle='None', alpha=0.6, label='individual cells'),
    mlines.Line2D([], [], color=CONTOUR_COLORS[0], linewidth=1.5,
                  label='densest 8% of area'),
    mlines.Line2D([], [], color=CONTOUR_COLORS[1], linewidth=2.5,
                  label='densest 1% of area'),
]
ax.legend(handles=handles, fontsize=9, loc='upper right', framealpha=0.92)

plt.tight_layout()
out_svg = f'{OUT_DIR}/Figure_4G_KDE_map.svg'
out_png = f'{OUT_DIR}/Figure_4G_KDE_map.png'
plt.savefig(out_svg, format='svg', bbox_inches='tight')
plt.savefig(out_png, dpi=300, bbox_inches='tight')
plt.close()

print(f'\nSaved: {out_svg}\n       {out_png}')
print(f'n = {n_cells:,} cells (report this in the figure legend)')
