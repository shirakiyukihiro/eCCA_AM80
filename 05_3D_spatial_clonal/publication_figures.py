#!/usr/bin/env python
# ============================================================
# 3D clonal spatial analysis — publication figures (KDE / Ripley / clone stats)
# Study      : Distinct CAF lineages and an Am80 stromal-reprogramming
#              strategy in extrahepatic cholangiocarcinoma (Meflin/ISLR).
# Repository : eCCA_AM80  (https://github.com/shirakiyukihiro/eCCA_AM80)
# Input = DoG clone-analysis .npz/.csv (see 3D README).
# Paths made relative for public release: put inputs in ./data, outputs in ./output.
# ============================================================

"""
Publication-quality figures for BDC tumor Meflin+ clonal analysis
====================================================================
Outputs 3 main figures as SVG (vector) + PNG (300 dpi):
  Figure 1: Detection methodology + uniformity validation
  Figure 2: Clonal expansion analysis (CORE)
  Figure 3: Per-clone characteristics

Styling: Nature/Cell-style, sans-serif, clean axes, journal column width
"""
import os
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.spatial import ConvexHull, cKDTree, Delaunay
from scipy.spatial.distance import pdist
from collections import Counter

# ============================================================
# Configuration
# ============================================================
NPZ_PATH  = 'data/BDC_tumor_DoG_Zcrop_clones.npz'
CSV_PATH  = 'data/BDC_tumor_DoG_Zcrop_cluster_features.csv'
OUT_DIR   = 'output'

# DBSCAN parameters (matched to original analysis)
EPS_UM        = 50
MIN_SAMPLES   = 3
N_MONTECARLO  = 100

# ============================================================
# Publication style
# ============================================================
plt.rcParams.update({
    'font.family':       'sans-serif',
    'font.sans-serif':   ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size':         8,
    'axes.labelsize':    9,
    'axes.titlesize':    10,
    'axes.labelweight':  'normal',
    'axes.titleweight':  'bold',
    'axes.spines.top':   False,
    'axes.spines.right': False,
    'axes.linewidth':    0.6,
    'xtick.direction':   'out',
    'ytick.direction':   'out',
    'xtick.major.size':  3,
    'ytick.major.size':  3,
    'xtick.major.width': 0.6,
    'ytick.major.width': 0.6,
    'xtick.labelsize':   7.5,
    'ytick.labelsize':   7.5,
    'legend.fontsize':   7.5,
    'legend.frameon':    False,
    'savefig.dpi':       300,
    'savefig.bbox':      'tight',
    'savefig.pad_inches': 0.05,
    'pdf.fonttype':      42,
    'svg.fonttype':      'none',
})

# Color palette
C_OBS    = '#1A535C'   # deep teal - observed/main data
C_NULL   = '#B8B8B8'   # neutral gray - null distribution
C_CLONE  = '#E63946'   # accent red - clones
C_GRAY   = '#9CA3AF'   # light gray - non-clone clusters
C_BLUE   = '#457B9D'   # blue - secondary highlights

PANEL_LABEL_KW = dict(fontsize=11, fontweight='bold', va='top', ha='left',
                       family='sans-serif')


def add_panel(ax, label, x=-0.15, y=1.05):
    # Axes3D requires text2D() for axes-fraction text placement
    if hasattr(ax, 'text2D'):
        ax.text2D(x, y, label, transform=ax.transAxes, **PANEL_LABEL_KW)
    else:
        ax.text(x, y, label, transform=ax.transAxes, **PANEL_LABEL_KW)


# ============================================================
# DBSCAN (fallback if sklearn missing)
# ============================================================
try:
    from sklearn.cluster import DBSCAN
    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False


def cluster_points(pts, eps, min_samples):
    if HAS_SKLEARN:
        return DBSCAN(eps=eps, min_samples=min_samples).fit(pts).labels_
    # Fallback
    n = len(pts); tree = cKDTree(pts)
    labels = -1 * np.ones(n, dtype=int)
    visited = np.zeros(n, dtype=bool); cid = 0
    for i in range(n):
        if visited[i]: continue
        visited[i] = True
        neighbors = tree.query_ball_point(pts[i], eps)
        if len(neighbors) < min_samples: continue
        labels[i] = cid
        seed = list(neighbors); idx = 0
        while idx < len(seed):
            j = seed[idx]; idx += 1
            if not visited[j]:
                visited[j] = True
                j_neigh = tree.query_ball_point(pts[j], eps)
                if len(j_neigh) >= min_samples:
                    for k in j_neigh:
                        if k not in seed: seed.append(k)
            if labels[j] == -1: labels[j] = cid
        cid += 1
    return labels


# ============================================================
# Load data
# ============================================================
print('Loading data...')
d = np.load(NPZ_PATH)
centroids = d['centroids_um'].astype(np.float64)  # (N, 3) z, y, x
labels = d['cluster_labels'].astype(int)
clone_thr = int(d['clone_threshold'][0])
n_total = len(centroids)
print(f'  Cells: {n_total}, clone threshold: ≥{clone_thr}')

# Features from CSV
features = []
with open(CSV_PATH) as f:
    for row in csv.DictReader(f):
        features.append({
            'cluster_id': int(row['cluster_id']),
            'n_cells':    int(row['n_cells']),
            'centroid_z': float(row['centroid_z']),
            'centroid_y': float(row['centroid_y']),
            'centroid_x': float(row['centroid_x']),
            'rms_radius_um': float(row['rms_radius_um']),
            'sphericity': float(row['sphericity']),
            'hull_volume_um3': float(row['hull_volume_um3']),
            'local_density_per_mm3': float(row['local_density_per_mm3']),
            'is_clone':  row['is_clone'] == 'True',
        })
features.sort(key=lambda x: -x['n_cells'])
clone_features = [f for f in features if f['is_clone']]
n_clusters = len(features)
n_clones = len(clone_features)
print(f'  Clusters: {n_clusters}, clones: {n_clones}')

sizes_arr = np.array([f['n_cells'] for f in features])

# Derived stats
hull = ConvexHull(centroids)
hull_vol_um3 = hull.volume
hull_vol_mm3 = hull_vol_um3 / 1e9

tree = cKDTree(centroids)
dd, _ = tree.query(centroids, k=2)
nn_dist = dd[:, 1]
mean_nn = nn_dist.mean()
expected_nn = 0.5 / (n_total / hull_vol_um3)**(1/3)
clark_evans = mean_nn / expected_nn
print(f'  Hull volume: {hull_vol_mm3:.2f} mm³, CE R: {clark_evans:.3f}')

# Inter-clone distance
clone_xyz = np.array([[f['centroid_z'], f['centroid_y'], f['centroid_x']]
                      for f in clone_features])
inter_d = pdist(clone_xyz)

# Non-clone mask for 3D plots
clone_ids = set(f['cluster_id'] for f in clone_features)
non_clone_mask = ~np.isin(labels, list(clone_ids))

# ============================================================
# Monte Carlo null
# ============================================================
print(f'Running Monte Carlo null ({N_MONTECARLO} iterations)...')
delaunay = Delaunay(centroids[hull.vertices])
bbox_min, bbox_max = centroids.min(axis=0), centroids.max(axis=0)
null_sizes_all = []
np.random.seed(42)
for sim in range(N_MONTECARLO):
    sampled = []
    while len(sampled) < n_total:
        cand = np.random.uniform(bbox_min, bbox_max, size=(n_total*3, 3))
        sampled.extend(cand[delaunay.find_simplex(cand) >= 0].tolist())
    rand_pts = np.array(sampled[:n_total])
    rl = cluster_points(rand_pts, EPS_UM, MIN_SAMPLES)
    valid = rl != -1
    if valid.sum() > 0:
        for c in Counter(rl[valid]).values():
            null_sizes_all.append(c)
    if (sim + 1) % 25 == 0:
        print(f'  {sim+1}/{N_MONTECARLO}')
null_sizes_all = np.array(null_sizes_all)


# ============================================================
# FIGURE 1: Detection & Validation
# ============================================================
print('\nFigure 1: Detection & Validation...')
fig1 = plt.figure(figsize=(7.2, 6.0))
gs1 = GridSpec(2, 2, figure=fig1, height_ratios=[1.0, 1.0],
                hspace=0.42, wspace=0.32, left=0.08, right=0.97,
                top=0.94, bottom=0.10)

# Panel A: 3D scatter of all cells, depth-coded
ax_a = fig1.add_subplot(gs1[0, 0], projection='3d')
sc = ax_a.scatter(centroids[:, 2], centroids[:, 1], centroids[:, 0],
                  c=centroids[:, 0], cmap='viridis',
                  s=0.6, alpha=0.6, edgecolors='none')
ax_a.set_xlabel('X (µm)', fontsize=7.5, labelpad=-3)
ax_a.set_ylabel('Y (µm)', fontsize=7.5, labelpad=-3)
ax_a.set_zlabel('Z (µm)', fontsize=7.5, labelpad=-3)
ax_a.tick_params(labelsize=6)
ax_a.view_init(elev=25, azim=-60)
ax_a.set_title(f'All detected cells\n(n = {n_total:,})', fontsize=9, fontweight='bold', pad=8)

# Panel B: Z-depth uniformity profile
ax_b = fig1.add_subplot(gs1[0, 1])
z_bins = np.linspace(0, 1500, 16)
z_counts, _ = np.histogram(centroids[:, 0], bins=z_bins)
zc = (z_bins[:-1] + z_bins[1:]) / 2
ax_b.bar(zc, z_counts, width=(z_bins[1] - z_bins[0]) * 0.88,
         color=C_OBS, edgecolor='white', linewidth=0.4)
ax_b.axhline(z_counts.mean(), color=C_CLONE, linestyle='--', linewidth=0.9,
             label=f'mean = {z_counts.mean():.0f}')
ax_b.set_xlabel('Tissue depth Z (µm)')
ax_b.set_ylabel('Cells per 100-µm Z bin')
ax_b.set_title('Detection uniformity', fontsize=9, fontweight='bold')
ax_b.set_xlim(-50, 1550)
ax_b.legend(loc='upper right')

# Panel C: Nearest-neighbor distance distribution
ax_c = fig1.add_subplot(gs1[1, 0])
bins_nn = np.linspace(0, 200, 41)
ax_c.hist(nn_dist, bins=bins_nn, color=C_OBS, edgecolor='white',
          linewidth=0.3, alpha=0.85)
ax_c.axvline(mean_nn, color=C_CLONE, linestyle='--', linewidth=0.9,
             label=f'observed mean: {mean_nn:.1f} µm')
ax_c.axvline(expected_nn, color=C_BLUE, linestyle=':', linewidth=1.2,
             label=f'CSR expected: {expected_nn:.1f} µm')
ax_c.set_xlabel('Nearest-neighbor distance (µm)')
ax_c.set_ylabel('Cell count')
ax_c.set_xlim(0, 200)
ax_c.set_title(f'Spatial clustering (Clark–Evans R = {clark_evans:.3f})',
               fontsize=9, fontweight='bold')
ax_c.legend(loc='upper right')

# Panel D: Summary statistics
ax_d = fig1.add_subplot(gs1[1, 1])
ax_d.axis('off')
n_in_clusters = sum(sizes_arr)
summary = (
    f'Detection summary\n'
    f'\n'
    f'Method                 DoG (3D blob detection)\n'
    f'σ_small / σ_large      5.3 / 10.6 µm\n'
    f'\n'
    f'Cells detected         {n_total:,}\n'
    f'Tumor volume           {hull_vol_mm3:.2f} mm³\n'
    f'Cell density           {n_total/hull_vol_mm3:.0f} cells/mm³\n'
    f'\n'
    f'Spatial test (CSR null)\n'
    f'Mean NN                {mean_nn:.1f} µm\n'
    f'Expected               {expected_nn:.1f} µm\n'
    f'Clark–Evans R          {clark_evans:.3f}\n'
    f'                       (R < 1: clustered)\n'
    f'\n'
    f'DBSCAN (ε = 50 µm)\n'
    f'Total clusters         {n_clusters}\n'
    f'Cells in clusters      {n_in_clusters:,}  ({100*n_in_clusters/n_total:.1f}%)\n'
    f'\n'
    f'Monte Carlo null (100 sim)\n'
    f'Clone threshold        ≥ {clone_thr} cells\n'
    f'                       (99th %ile)\n'
)
ax_d.text(0.0, 0.96, summary, transform=ax_d.transAxes,
          fontsize=7.5, family='monospace', va='top', linespacing=1.5)

add_panel(ax_a, 'A', x=-0.05, y=1.02)
add_panel(ax_b, 'B')
add_panel(ax_c, 'C')
add_panel(ax_d, 'D', x=0.0)

fig1.savefig(f'{OUT_DIR}/Figure_1_detection_validation.svg', format='svg')
fig1.savefig(f'{OUT_DIR}/Figure_1_detection_validation.png', dpi=300)
plt.close(fig1)
print('  Saved Figure 1')


# ============================================================
# FIGURE 2: Clonal analysis (CORE)
# ============================================================
print('\nFigure 2: Clonal analysis...')
fig2 = plt.figure(figsize=(7.2, 7.5))
gs2 = GridSpec(3, 3, figure=fig2, height_ratios=[1.3, 1.0, 0.9],
                hspace=0.5, wspace=0.42, left=0.07, right=0.96,
                top=0.96, bottom=0.05)

# Panel A: 3D scatter, top 20 clones colored
ax_a = fig2.add_subplot(gs2[0, 0:2], projection='3d')
ax_a.scatter(centroids[non_clone_mask, 2], centroids[non_clone_mask, 1],
             centroids[non_clone_mask, 0],
             s=0.5, alpha=0.15, c='lightgray', edgecolors='none')
clone_cmap = plt.cm.tab20
for i, cf in enumerate(clone_features[:20]):
    m = labels == cf['cluster_id']
    ax_a.scatter(centroids[m, 2], centroids[m, 1], centroids[m, 0],
                 s=12, alpha=0.85, color=clone_cmap(i % 20),
                 edgecolor='white', linewidth=0.25)
ax_a.set_xlabel('X (µm)', fontsize=7.5, labelpad=-2)
ax_a.set_ylabel('Y (µm)', fontsize=7.5, labelpad=-2)
ax_a.set_zlabel('Z (µm)', fontsize=7.5, labelpad=-2)
ax_a.tick_params(labelsize=6)
ax_a.view_init(elev=22, azim=-55)
ax_a.set_title(f'Spatial distribution of top 20 clones\n'
               f'(of {n_clones} total clones, gray = non-clone)',
               fontsize=9, fontweight='bold', pad=8)

# Panel B: Cluster size: observed vs null
ax_b = fig2.add_subplot(gs2[0, 2])
max_size = max(sizes_arr.max(), null_sizes_all.max())
bins = np.logspace(np.log10(2.5), np.log10(max_size + 1), 22)
bc = (bins[:-1] + bins[1:]) / 2
bw = np.diff(bins) * 0.45
null_per_sim, _ = np.histogram(null_sizes_all, bins=bins)
null_per_sim = null_per_sim / N_MONTECARLO
obs_counts, _ = np.histogram(sizes_arr, bins=bins)
ax_b.bar(bc, obs_counts, width=np.diff(bins) * 0.85,
         color=C_OBS, alpha=0.95, label='Observed',
         edgecolor='white', linewidth=0.3, zorder=2)
ax_b.bar(bc, null_per_sim, width=np.diff(bins) * 0.85,
         color=C_NULL, alpha=0.75, label='CSR null (mean)',
         edgecolor='white', linewidth=0.3, zorder=1)
ax_b.axvline(clone_thr, color=C_CLONE, linestyle='--', linewidth=0.9,
             label=f'Clone threshold = {clone_thr}')
ax_b.set_xscale('log'); ax_b.set_yscale('log')
ax_b.set_xlabel('Cluster size (cells)')
ax_b.set_ylabel('Cluster count')
ax_b.set_title('Observed vs random null', fontsize=9, fontweight='bold')
ax_b.legend(loc='upper right')

# Panel C: Top 20 clones bar chart
ax_c = fig2.add_subplot(gs2[1, 0:2])
top20 = sizes_arr[:20]
colors = [C_CLONE if s >= clone_thr else C_GRAY for s in top20]
ax_c.bar(range(1, 21), top20, color=colors, edgecolor='white', linewidth=0.5)
ax_c.axhline(clone_thr, color=C_BLUE, linestyle='--', linewidth=0.9,
             label=f'Clone threshold = {clone_thr}')
ax_c.set_xticks(range(1, 21))
ax_c.set_xlabel('Cluster rank')
ax_c.set_ylabel('Cluster size (cells)')
ax_c.set_title(f'Largest 20 clusters (out of {n_clusters} total)',
               fontsize=9, fontweight='bold')
ax_c.legend(loc='upper right')
for i, v in enumerate(top20[:10]):
    ax_c.text(i + 1, v + 1.5, str(int(v)), ha='center', va='bottom', fontsize=6.5)

# Panel D: Inter-clone distance
ax_d = fig2.add_subplot(gs2[1, 2])
ax_d.hist(inter_d / 1000, bins=22, color=C_OBS, edgecolor='white',
          linewidth=0.3, alpha=0.9)
ax_d.axvline(np.median(inter_d) / 1000, color=C_CLONE, linestyle='--',
             linewidth=0.9, label=f'median: {np.median(inter_d)/1000:.2f} mm')
ax_d.set_xlabel('Inter-clone distance (mm)')
ax_d.set_ylabel('Pair count')
ax_d.set_title(f'Multifocality\n({len(inter_d):,} pairs from {n_clones} clones)',
               fontsize=9, fontweight='bold')
ax_d.legend(loc='upper right')

# Panel E: Top 10 clones table
ax_e = fig2.add_subplot(gs2[2, :])
ax_e.axis('off')
table_data = [['Rank', 'Cluster\nID', 'Cells', 'Z (µm)', 'RMS radius\n(µm)',
                'Sphericity', 'Local density\n(cells/mm³)']]
for i, f in enumerate(clone_features[:10]):
    table_data.append([
        f'{i+1}', f"{f['cluster_id']}", f"{f['n_cells']}",
        f"{f['centroid_z']:.0f}", f"{f['rms_radius_um']:.0f}",
        f"{f['sphericity']:.3f}", f"{f['local_density_per_mm3']:,.0f}",
    ])
table = ax_e.table(cellText=table_data, loc='center', cellLoc='center',
                    colWidths=[0.07, 0.10, 0.08, 0.10, 0.13, 0.12, 0.16])
table.auto_set_font_size(False)
table.set_fontsize(7.5)
table.scale(1, 1.45)
for j in range(len(table_data[0])):
    table[(0, j)].set_facecolor('#E8EEF1')
    table[(0, j)].set_text_props(weight='bold')
# Highlight largest clone
for j in range(len(table_data[0])):
    table[(1, j)].set_facecolor('#FFF3E6')
ax_e.set_title('Top 10 clones — quantitative characteristics',
               fontsize=9, fontweight='bold', pad=4)

add_panel(ax_a, 'A', x=-0.02, y=1.0)
add_panel(ax_b, 'B', x=-0.28)
add_panel(ax_c, 'C', x=-0.06)
add_panel(ax_d, 'D', x=-0.28)
add_panel(ax_e, 'E', x=-0.02, y=1.04)

fig2.savefig(f'{OUT_DIR}/Figure_2_clonal_analysis.svg', format='svg')
fig2.savefig(f'{OUT_DIR}/Figure_2_clonal_analysis.png', dpi=300)
plt.close(fig2)
print('  Saved Figure 2')


# ============================================================
# FIGURE 3: Clone characteristics
# ============================================================
print('\nFigure 3: Clone characteristics...')
fig3 = plt.figure(figsize=(7.2, 6.2))
gs3 = GridSpec(2, 2, figure=fig3, hspace=0.45, wspace=0.35,
                left=0.09, right=0.97, top=0.94, bottom=0.10)

# Panel A: Sphericity vs size
ax_a = fig3.add_subplot(gs3[0, 0])
sizes_all = [f['n_cells'] for f in features]
sph_all = [f['sphericity'] for f in features]
cmask = [f['is_clone'] for f in features]
nc_x = [s for s, c in zip(sizes_all, cmask) if not c]
nc_y = [p for p, c in zip(sph_all, cmask) if not c]
c_x  = [s for s, c in zip(sizes_all, cmask) if c]
c_y  = [p for p, c in zip(sph_all, cmask) if c]
ax_a.scatter(nc_x, nc_y, s=10, alpha=0.4, c=C_GRAY,
             edgecolors='none', label='Non-clone clusters')
ax_a.scatter(c_x, c_y, s=22, alpha=0.85, c=C_CLONE,
             edgecolor='white', linewidth=0.4, label=f'Clones (n = {n_clones})')
ax_a.set_xscale('log')
ax_a.set_xlabel('Cluster size (cells)')
ax_a.set_ylabel('Sphericity (0 = line, 1 = sphere)')
ax_a.set_ylim(-0.02, 0.45)
ax_a.set_title('Morphology vs size', fontsize=9, fontweight='bold')
ax_a.legend(loc='upper left')

# Panel B: Density vs size (n>=5 only to avoid hull artifacts)
ax_b = fig3.add_subplot(gs3[0, 1])
filt = [(f['n_cells'], f['local_density_per_mm3'], f['is_clone'])
        for f in features if f['n_cells'] >= 5 and f['local_density_per_mm3'] > 0]
if filt:
    fs, fd, fc = zip(*filt)
    nc_x = [s for s, c in zip(fs, fc) if not c]
    nc_y = [d for d, c in zip(fd, fc) if not c]
    c_x  = [s for s, c in zip(fs, fc) if c]
    c_y  = [d for d, c in zip(fd, fc) if c]
    ax_b.scatter(nc_x, nc_y, s=10, alpha=0.4, c=C_GRAY,
                 edgecolors='none', label='Non-clone')
    ax_b.scatter(c_x, c_y, s=22, alpha=0.85, c=C_CLONE,
                 edgecolor='white', linewidth=0.4, label='Clones')
ax_b.set_xscale('log'); ax_b.set_yscale('log')
ax_b.set_xlabel('Cluster size (cells)')
ax_b.set_ylabel('Local density (cells/mm³)')
ax_b.set_title('Density vs size (clusters n ≥ 5)', fontsize=9, fontweight='bold')
ax_b.legend(loc='lower left')

# Panel C: 3D distribution of clone centroids
ax_c = fig3.add_subplot(gs3[1, 0], projection='3d')
ax_c.scatter(centroids[non_clone_mask, 2], centroids[non_clone_mask, 1],
             centroids[non_clone_mask, 0],
             s=0.3, alpha=0.07, c='lightgray', edgecolors='none')
clone_xs = [f['centroid_x'] for f in clone_features]
clone_ys = [f['centroid_y'] for f in clone_features]
clone_zs = [f['centroid_z'] for f in clone_features]
clone_n  = np.array([f['n_cells'] for f in clone_features])
sc3 = ax_c.scatter(clone_xs, clone_ys, clone_zs,
                    s=clone_n * 2.5, c=clone_n, cmap='magma',
                    alpha=0.85, edgecolor='white', linewidth=0.3,
                    vmin=clone_thr, vmax=75)
cb = plt.colorbar(sc3, ax=ax_c, shrink=0.55, pad=0.10, aspect=14)
cb.set_label('Cells per clone', fontsize=7)
cb.ax.tick_params(labelsize=6)
ax_c.set_xlabel('X (µm)', fontsize=7.5, labelpad=-2)
ax_c.set_ylabel('Y (µm)', fontsize=7.5, labelpad=-2)
ax_c.set_zlabel('Z (µm)', fontsize=7.5, labelpad=-2)
ax_c.tick_params(labelsize=6)
ax_c.view_init(elev=20, azim=-50)
ax_c.set_title(f'Spatial distribution of {n_clones} clone foci',
                fontsize=9, fontweight='bold', pad=8)

# Panel D: Reverse cumulative distribution
ax_d = fig3.add_subplot(gs3[1, 1])
sorted_sizes = np.sort(sizes_arr)
cdf_rev = np.arange(len(sorted_sizes), 0, -1)
ax_d.loglog(sorted_sizes, cdf_rev, 'o-', color=C_OBS, markersize=3.0,
            markerfacecolor=C_OBS, markeredgecolor='white',
            markeredgewidth=0.4, linewidth=0.8)
ax_d.axvline(clone_thr, color=C_CLONE, linestyle='--', linewidth=0.9,
             label=f'Clone threshold = {clone_thr}')
ax_d.set_xlabel('Cluster size ≥ N')
ax_d.set_ylabel('Number of clusters')
ax_d.set_title('Reverse cumulative distribution', fontsize=9, fontweight='bold')
ax_d.legend(loc='upper right')

add_panel(ax_a, 'A')
add_panel(ax_b, 'B')
add_panel(ax_c, 'C', x=-0.04, y=1.0)
add_panel(ax_d, 'D')

fig3.savefig(f'{OUT_DIR}/Figure_3_clone_characteristics.svg', format='svg')
fig3.savefig(f'{OUT_DIR}/Figure_3_clone_characteristics.png', dpi=300)
plt.close(fig3)
print('  Saved Figure 3')


# ============================================================
# Summary
# ============================================================
print()
print('=' * 70)
print('  PUBLICATION FIGURES GENERATED')
print('=' * 70)
print()
print('Output files (vector SVG + raster PNG @ 300 dpi):')
print(f'  {OUT_DIR}/Figure_1_detection_validation.svg/.png')
print(f'  {OUT_DIR}/Figure_2_clonal_analysis.svg/.png')
print(f'  {OUT_DIR}/Figure_3_clone_characteristics.svg/.png')
print()
print('For journal submission: use SVG files.')
print('SVG text is preserved as text → editable in Illustrator/Inkscape.')
print('PNG files are for preview/review.')
print('=' * 70)
