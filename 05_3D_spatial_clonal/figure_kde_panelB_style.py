#!/usr/bin/env python
# ============================================================
# Figure 4G — spatial map of tdTomato+ Meflin-lineage cell clusters (KDE panel)
# Study      : Distinct CAF lineages and an Am80 stromal-reprogramming strategy
#              in extrahepatic cholangiocarcinoma (Meflin/ISLR).
# Repository : eCCA_AM80  (https://github.com/shirakiyukihiro/eCCA_AM80)
# Purpose    : Renders the published Fig. 4G panel: every detected cell as a point,
#              two-level kernel-density contours (92nd / 99th percentile) marking
#              progressively higher-density regions, and the ten largest DBSCAN
#              clusters (eps = 50 um) numbered at their centroids.
# Input      : output/BDC_tumor_DoG_Zcrop_clones.npz  (from analyze_tumor_clones_DoG_Zcrop.py)
# Output     : output/Figure_KDE_panelB_style.svg / .png
# ============================================================

"""
過去 (FINAL v10) Panel B style を現在 (DoG Z<1500) データで再現
==================================================================

Panel B 構成:
  - 灰色 scatter (全細胞)
  - KDE 2 段階 contour line (92%, 99% percentile)
  - Top 10 clones: 黄色塗り + 赤縁 marker + 黒番号
  - 統合 inset (cluster region | focal cluster | noise の概念図)
  - 凡例
"""
import os, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.ndimage import label as nd_label
import matplotlib.lines as mlines
from matplotlib.patches import Ellipse, Polygon
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from collections import Counter

# ============================================================
# Configuration
# ============================================================
NPZ_PATH = 'output/BDC_tumor_DoG_Zcrop_clones.npz'
OUT_DIR  = 'output'

KDE_BW            = 0.03   # KDE bandwidth (smaller = sharper contours)
CONTOUR_PCTS      = [92, 99]  # contour percentiles
TOP_N_HIGHLIGHT   = 10

# ============================================================
# Load data
# ============================================================
print('Loading data...')
d = np.load(NPZ_PATH)
cents = d['centroids_um'].astype(np.float64)
labels = d['cluster_labels'].astype(int)
N = len(cents)

cl_u, cl_c = np.unique(labels[labels != -1], return_counts=True)
order = np.argsort(-cl_c)
cl_u_sorted = cl_u[order]
cl_c_sorted = cl_c[order]
top_ids = cl_u_sorted[:TOP_N_HIGHLIGHT]
top_sizes = cl_c_sorted[:TOP_N_HIGHLIGHT]
top_centers = np.array([cents[labels == cid].mean(axis=0) for cid in top_ids])

clone_thr = int(d['clone_threshold'][0])
n_clones = sum(1 for c in cl_c if c >= clone_thr)
print(f'  N={N}, clusters={len(cl_u)}, clones (>={clone_thr})={n_clones}')
print(f'  Top {TOP_N_HIGHLIGHT} sizes: {top_sizes.tolist()}')

# ============================================================
# Figure: Panel B style (standalone)
# ============================================================
print('Generating figure...')
fig, ax_B = plt.subplots(figsize=(8, 8))

# Background: all cells in gray (light)
ax_B.scatter(cents[:, 2], cents[:, 1], s=3, color='gray', alpha=0.35, zorder=1)

# KDE compute on XY projection
xy = cents[:, [2, 1]]
kde = gaussian_kde(xy.T, bw_method=KDE_BW)
x_range = (xy[:, 0].min() - 50, xy[:, 0].max() + 50)
y_range = (xy[:, 1].min() - 50, xy[:, 1].max() + 50)
X, Y = np.meshgrid(np.linspace(*x_range, 250), np.linspace(*y_range, 250))
Z = kde(np.vstack([X.ravel(), Y.ravel()])).reshape(X.shape)

# 2-level contour lines (NOT heatmap)
levels = [np.percentile(Z, p) for p in CONTOUR_PCTS]
contour_colors_c = ['#FFB6C1', '#DC143C']  # light pink, crimson
contour_widths = [1.2, 2.2]
for lv, cc, lw in zip(levels, contour_colors_c, contour_widths):
    try:
        ax_B.contour(X, Y, Z, levels=[lv], colors=[cc], linewidths=lw)
    except Exception:
        pass

# Count high-density spots (>95 percentile)
n_spots = nd_label(Z > np.percentile(Z, 95))[1]

# Top 10 cluster centers: yellow filled + red edge + black number
for i, (center, sz) in enumerate(zip(top_centers, top_sizes)):
    ax_B.scatter(center[2], center[1], s=160,
                  edgecolor='red', facecolor='yellow',
                  linewidth=2, alpha=0.85, zorder=10)
    ax_B.text(center[2], center[1], f'{i + 1}',
              color='black', fontsize=9,
              ha='center', va='center', fontweight='bold', zorder=11)

# Axes
ax_B.set_xlabel('X (µm)', fontsize=12)
ax_B.set_ylabel('Y (µm)', fontsize=12)
ax_B.set_aspect('equal')
ax_B.invert_yaxis()
ax_B.grid(alpha=0.3)
ax_B.set_title('Spatial map of tdTomato⁺ cell clusters across the tumor\n'
                f'(n = {N:,} cells, {n_clones} clones)',
               fontsize=13, fontweight='bold', loc='left', pad=10)

# ============================================================
# Integrated INSET: 3-state conceptual diagram
# ============================================================
inset_ax = inset_axes(ax_B, width='55%', height='34%',
                       loc='lower left', borderpad=1.5)
inset_ax.set_xlim(0, 1); inset_ax.set_ylim(0, 1)
from matplotlib.patches import Rectangle
inset_ax.add_patch(Rectangle((0, 0), 1, 1, facecolor='white',
                                edgecolor='none', zorder=0))
for spine in inset_ax.spines.values():
    spine.set_edgecolor('black')
    spine.set_linewidth(1.2)
inset_ax.set_xticks([]); inset_ax.set_yticks([])
inset_ax.set_zorder(100)

# LEFT: "cluster region" - loose cells with light pink outline
loose_cells = [(0.07, 0.55), (0.15, 0.62), (0.20, 0.50),
                (0.27, 0.60), (0.32, 0.48)]
loose_outline = [(0.04, 0.55), (0.10, 0.70), (0.20, 0.72), (0.30, 0.68),
                  (0.35, 0.58), (0.32, 0.42), (0.20, 0.40), (0.08, 0.45)]
inset_ax.add_patch(Polygon(loose_outline, closed=True, facecolor='none',
                              edgecolor='#FFB6C1', linewidth=2.5, zorder=4))
loose_pairs = [(0, 1), (0, 2), (1, 3), (2, 3), (3, 4), (2, 4)]
for i, j in loose_pairs:
    inset_ax.plot([loose_cells[i][0], loose_cells[j][0]],
                    [loose_cells[i][1], loose_cells[j][1]],
                    '-', color='lightgray', lw=0.8, zorder=3)
for cx, cy in loose_cells:
    inset_ax.scatter(cx, cy, s=70, color='red',
                       edgecolor='black', linewidth=1, zorder=10)
inset_ax.text(0.20, 0.88, 'cluster region',
                fontsize=10, ha='center', fontweight='bold', color='#C71585')
inset_ax.text(0.20, 0.22, 'cells loosely\nlinked (≤50 µm)',
                fontsize=8, ha='center', style='italic', color='dimgray')

# separator
inset_ax.plot([0.42, 0.42], [0.10, 0.94], color='gray',
                linestyle=':', lw=1.3, alpha=0.6)

# MIDDLE: "focal cluster" - densely packed cells with crimson outline
focal_cells = [(0.56, 0.56), (0.60, 0.62), (0.64, 0.55), (0.59, 0.49),
                (0.65, 0.62), (0.68, 0.55), (0.62, 0.58)]
focal_outline = [(0.51, 0.55), (0.55, 0.66), (0.65, 0.68), (0.72, 0.60),
                  (0.72, 0.50), (0.65, 0.44), (0.55, 0.46)]
inset_ax.add_patch(Polygon(focal_outline, closed=True, facecolor='none',
                              edgecolor='#DC143C', linewidth=3, zorder=4))
for cx, cy in focal_cells:
    inset_ax.scatter(cx, cy, s=70, color='red',
                       edgecolor='black', linewidth=1, zorder=10)
inset_ax.text(0.62, 0.88, 'focal cluster',
                fontsize=10, ha='center', fontweight='bold', color='#DC143C')
inset_ax.text(0.62, 0.22, 'cells densely\npacked',
                fontsize=8, ha='center', style='italic', color='dimgray')

# separator
inset_ax.plot([0.78, 0.78], [0.10, 0.94], color='gray',
                linestyle=':', lw=1.3, alpha=0.6)

# RIGHT: "noise" - single isolated cell
fig_w, fig_h = fig.get_size_inches()
ax_pos = inset_ax.get_position()
inset_aspect = (ax_pos.height * fig_h) / (ax_pos.width * fig_w)
rx = 0.025; ry = rx / inset_aspect

inset_ax.add_patch(Ellipse((0.89, 0.55), 2*rx, 2*ry,
                              facecolor='lightgray', alpha=0.3,
                              edgecolor='gray', linewidth=0.7, zorder=2))
inset_ax.scatter(0.89, 0.55, s=70, color='gray',
                   edgecolor='black', linewidth=1, zorder=10)
inset_ax.text(0.89, 0.88, 'noise',
                fontsize=10, ha='center', fontweight='bold', color='dimgray')
inset_ax.text(0.89, 0.22, 'isolated\ncell',
                fontsize=8, ha='center', style='italic', color='dimgray')

# Scale bar
inset_ax.annotate('', xy=(0.15, 0.08), xytext=(0.05, 0.08),
                    arrowprops=dict(arrowstyle='<->', color='black', lw=1.2))
inset_ax.text(0.10, 0.02, '50 µm', fontsize=8, ha='center', fontweight='bold')

# ============================================================
# Legend (top right)
# ============================================================
cell_marker = mlines.Line2D([], [], color='gray', marker='o', markersize=6,
                              linestyle='None', alpha=0.6,
                              label='individual cells')
contour_92_handle = mlines.Line2D([], [], color='#FFB6C1', linewidth=1.5,
                                    label='92nd %ile contour')
contour_99_handle = mlines.Line2D([], [], color='#DC143C', linewidth=2.5,
                                    label='99th %ile contour')
center_marker = mlines.Line2D([], [], color='red', marker='o', markersize=10,
                                markerfacecolor='yellow', linestyle='None',
                                label='Top 10 (numbered)')
ax_B.legend(handles=[cell_marker, contour_92_handle, contour_99_handle, center_marker],
             fontsize=9, loc='upper right', framealpha=0.92)

# ============================================================
# Save
# ============================================================
plt.tight_layout()
out_svg = f'{OUT_DIR}/Figure_KDE_panelB_style.svg'
out_png = f'{OUT_DIR}/Figure_KDE_panelB_style.png'
plt.savefig(out_svg, format='svg', bbox_inches='tight')
plt.savefig(out_png, dpi=300, bbox_inches='tight')
plt.close()

print()
print('=' * 60)
print(f'  Saved: {out_svg}')
print(f'         {out_png}')
print('=' * 60)
print()
print(f'Top {TOP_N_HIGHLIGHT} clones (numbered in figure):')
for i, (cid, sz, ctr) in enumerate(zip(top_ids, top_sizes, top_centers)):
    print(f'  #{i+1:2d}: id={cid:3d}, size={sz:3d}, '
            f'(X, Y) = ({ctr[2]:.0f}, {ctr[1]:.0f}) µm')
print()
print(f'High-density spots (>95% KDE): ~{n_spots}')
