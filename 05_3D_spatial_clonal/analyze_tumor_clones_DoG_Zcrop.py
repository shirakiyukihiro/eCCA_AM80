#!/usr/bin/env python
# ============================================================
# Spatial statistics and clustering (Clark-Evans + DBSCAN + Monte-Carlo null), Z < 1500 um
# Study      : Distinct CAF lineages and an Am80 stromal-reprogramming strategy
#              in extrahepatic cholangiocarcinoma (Meflin/ISLR).
# Repository : eCCA_AM80  (https://github.com/shirakiyukihiro/eCCA_AM80)
# Clusters detected cells (DBSCAN, eps=50 um), tests departure from complete spatial
# randomness, and derives a size threshold from a Monte-Carlo CSR null.
# Paths made relative for public release: raw .ims files in ./data (not distributed;
# available on request), region masks in ./data/mask_polygons, results in ./output.
# ============================================================

"""
Spatial statistics and clustering of detected cells (Z < 1500 um).

Light-sheet signal attenuates with depth, so the analysis is restricted to the
uniformly detected zone (Z < 1500 um; see validate_detection_uniformity.py).

Computes:
  - Clark-Evans index R (3D), with a Monte-Carlo p-value
  - DBSCAN clusters (eps = 50 um, min_samples = 3)
  - a cluster-size threshold from a Monte-Carlo complete-spatial-randomness null

The numbers reported in the Results text come from this script.
"""
import os, sys, json
from math import gamma, pi
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.spatial import cKDTree, ConvexHull, Delaunay
from collections import Counter

INPUT_NPZ      = 'output/result_BDC_tumor_DoG.npz'
OUTPUT_NPZ     = 'output/BDC_tumor_DoG_Zcrop_clones.npz'
OUTPUT_CSV     = 'output/BDC_tumor_DoG_Zcrop_cluster_features.csv'
OUTPUT_PNG     = 'output/BDC_tumor_DoG_Zcrop_clone_analysis.png'
OUTPUT_SVG     = 'output/BDC_tumor_DoG_Zcrop_clone_analysis.svg'

Z_MAX_UM       = 1500   # ★ Z 上限 (uniform detection 範囲)
EPS_UM         = 50
MIN_SAMPLES    = 3
N_MONTECARLO   = 100
CLONE_PCTILE   = 99

try:
    from sklearn.cluster import DBSCAN
    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False


def dbscan_kdtree(pts, eps, min_samples):
    n = len(pts)
    tree = cKDTree(pts)
    labels = -1 * np.ones(n, dtype=int)
    visited = np.zeros(n, dtype=bool)
    cid = 0
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


def cluster_points(pts, eps, min_samples):
    if HAS_SKLEARN:
        return DBSCAN(eps=eps, min_samples=min_samples).fit(pts).labels_
    return dbscan_kdtree(pts, eps, min_samples)


print('=' * 70)
print(f'  BDC tumor (DoG, Z < {Z_MAX_UM}µm) clone analysis')
print('=' * 70)

# Load and filter
d = np.load(INPUT_NPZ)
centroids_all = d['centroids_um'].astype(np.float64)
n_before = len(centroids_all)
keep = centroids_all[:, 0] < Z_MAX_UM
centroids = centroids_all[keep]
n_total = len(centroids)
print(f'\nFiltered to Z < {Z_MAX_UM}µm: {n_total}/{n_before} cells ({100*n_total/n_before:.1f}%)')

# Stats
hull = ConvexHull(centroids)
hull_vol_um3 = hull.volume
hull_vol_mm3 = hull_vol_um3 / 1e9
print(f'3D hull volume: {hull_vol_mm3:.3f} mm³')
print(f'Density: {n_total/hull_vol_mm3:.1f} cells/mm³')

# DBSCAN
print(f'\nDBSCAN (eps={EPS_UM}µm, min={MIN_SAMPLES})...')
labels = cluster_points(centroids, EPS_UM, MIN_SAMPLES)
n_clusters = labels.max() + 1 if labels.max() >= 0 else 0
n_noise = (labels == -1).sum()
cluster_sizes = Counter(labels[labels != -1])
sizes_arr = np.array(sorted(cluster_sizes.values(), reverse=True))
print(f'  Clusters: {n_clusters}, in-cluster: {n_total-n_noise} ({100*(n_total-n_noise)/n_total:.1f}%)')
print(f'  Top 10: {sizes_arr[:10].tolist()}')

# Spatial pattern
tree = cKDTree(centroids)
dd, _ = tree.query(centroids, k=2)
nn_dist = dd[:, 1]
mean_nn = nn_dist.mean()
# Expected nearest-neighbour distance under CSR in 3D:
#   E[d] = Gamma(4/3) / ((4/3) * pi * rho)**(1/3) = 0.553961 * rho**(-1/3)
# (0.5 is the 2D coefficient and must not be used here.)
C_3D = gamma(4 / 3) / ((4 / 3) * pi) ** (1 / 3)          # = 0.553961
rho = n_total / hull_vol_um3                             # cells per um^3
expected_nn = C_3D * rho ** (-1 / 3)
clark_evans = mean_nn / expected_nn
print(f'\nMean NN: {mean_nn:.1f}µm, Expected CSR: {expected_nn:.1f}µm, CE R: {clark_evans:.3f}')

# Monte Carlo null
print(f'\nMonte Carlo null ({N_MONTECARLO} iterations)...')
delaunay = Delaunay(centroids[hull.vertices])
bbox_min = centroids.min(axis=0); bbox_max = centroids.max(axis=0)

null_max_sizes = []; null_n_clusters_arr = []; null_sizes_all = []
null_mean_nn = []
np.random.seed(42)
for sim in range(N_MONTECARLO):
    sampled = []
    while len(sampled) < n_total:
        cand = np.random.uniform(bbox_min, bbox_max, size=(n_total * 3, 3))
        inside = delaunay.find_simplex(cand) >= 0
        sampled.extend(cand[inside].tolist())
    rand_pts = np.array(sampled[:n_total])

    # null mean nearest-neighbour distance (edge effects handled implicitly,
    # because the points are drawn inside the same convex hull)
    rdd, _ = cKDTree(rand_pts).query(rand_pts, k=2)
    null_mean_nn.append(rdd[:, 1].mean())

    rl = cluster_points(rand_pts, EPS_UM, MIN_SAMPLES)
    valid = rl != -1
    if valid.sum() > 0:
        rc = Counter(rl[valid])
        null_max_sizes.append(max(rc.values()))
        null_n_clusters_arr.append(len(rc))
        null_sizes_all.extend(rc.values())
    if (sim + 1) % 10 == 0: print(f'  {sim+1}/{N_MONTECARLO}')

null_mean_nn = np.array(null_mean_nn)
R_mc = mean_nn / null_mean_nn.mean()                       # edge-corrected R
p_mc = (int((null_mean_nn <= mean_nn).sum()) + 1) / (N_MONTECARLO + 1)
print(f'\nMonte-Carlo CSR null ({N_MONTECARLO} iterations):')
print(f'  null mean NN       : {null_mean_nn.mean():.1f} +/- {null_mean_nn.std():.1f} um')
print(f'  observed mean NN   : {mean_nn:.1f} um')
print(f'  R (edge-corrected) : {R_mc:.3f}')
print(f'  p (one-sided)      : {p_mc:.3f}'
      + ('   -> report as p < 0.01' if p_mc <= 0.01 else ''))

null_sizes_all = np.array(null_sizes_all)
clone_threshold = int(np.percentile(null_sizes_all, CLONE_PCTILE))
print(f'\nNull {CLONE_PCTILE}%ile: {clone_threshold}')

clone_mask = sizes_arr >= clone_threshold
n_clones = int(clone_mask.sum())
print(f'★ CLONE CANDIDATES (size ≥ {clone_threshold}): {n_clones}')

# Features
features = []
for cid in range(n_clusters):
    mask_c = labels == cid
    pts = centroids[mask_c]
    n = len(pts)
    if n < 3: continue
    cent = pts.mean(axis=0)
    rms_r = np.sqrt(np.mean(np.linalg.norm(pts - cent, axis=1)**2))
    cov = np.cov(pts.T)
    eigvals = np.sort(np.linalg.eigvalsh(cov))
    sph = eigvals[0]/eigvals[-1] if eigvals[-1] > 1e-6 else 0
    hv = 0; dens = 0
    if n >= 5:
        try:
            ch = ConvexHull(pts)
            hv = ch.volume
            if hv > 0: dens = n / (hv * 1e-9)
        except Exception: pass
    features.append({
        'cluster_id': cid, 'n_cells': n,
        'centroid_z': cent[0], 'centroid_y': cent[1], 'centroid_x': cent[2],
        'rms_radius_um': rms_r, 'sphericity': sph,
        'hull_volume_um3': hv, 'local_density_per_mm3': dens,
        'is_clone': bool(n >= clone_threshold),
    })
features.sort(key=lambda x: -x['n_cells'])

# Inter-clone
clone_features = [f for f in features if f['is_clone']]
inter_d = np.array([])
if len(clone_features) >= 2:
    from scipy.spatial.distance import pdist
    cc = np.array([[f['centroid_z'], f['centroid_y'], f['centroid_x']] for f in clone_features])
    inter_d = pdist(cc)

# Save CSV
import csv
with open(OUTPUT_CSV, 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(features[0].keys())
    for r in features: w.writerow(r.values())

np.savez(OUTPUT_NPZ,
         centroids_um=centroids.astype(np.float32),
         cluster_labels=labels.astype(np.int32),
         clone_threshold=np.array([clone_threshold]))

# Plot (simplified)
fig = plt.figure(figsize=(15, 10))
gs = GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.35)

# Cluster sizes
ax = fig.add_subplot(gs[0, 0])
max_size = max(sizes_arr.max(), null_sizes_all.max())
bins = np.logspace(np.log10(3), np.log10(max_size + 1), 30)
ax.hist(sizes_arr, bins=bins, alpha=0.7, label='Observed', color='red', edgecolor='black')
ax.hist(null_sizes_all, bins=bins, alpha=0.4, label=f'Null', color='gray', edgecolor='black',
        weights=np.ones_like(null_sizes_all)/N_MONTECARLO)
ax.axvline(clone_threshold, color='blue', linestyle='--', label=f'Clone thr = {clone_threshold}')
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlabel('Cluster size'); ax.set_ylabel('Count')
ax.set_title('Cluster sizes (vs null)', fontweight='bold')
ax.legend(fontsize=8); ax.grid(alpha=0.3)

# Top 20
ax = fig.add_subplot(gs[0, 1])
top20 = sizes_arr[:20]
colors_b = ['red' if s >= clone_threshold else 'gray' for s in top20]
ax.bar(range(1, len(top20)+1), top20, color=colors_b, edgecolor='black')
ax.axhline(clone_threshold, color='blue', linestyle='--')
ax.set_xlabel('Cluster rank'); ax.set_ylabel('Size')
ax.set_title(f'Top 20 (red=clone)', fontweight='bold')
ax.grid(alpha=0.3, axis='y')

# Sphericity vs size
ax = fig.add_subplot(gs[0, 2])
sf = [f['n_cells'] for f in features]; pf = [f['sphericity'] for f in features]
cf = [f['is_clone'] for f in features]
ax.scatter([s for s,c in zip(sf,cf) if not c], [p for p,c in zip(pf,cf) if not c],
           s=15, alpha=0.3, c='gray', label='Non-clone')
ax.scatter([s for s,c in zip(sf,cf) if c], [p for p,c in zip(pf,cf) if c],
           s=30, alpha=0.8, c='red', edgecolor='black', label='Clone')
ax.set_xscale('log'); ax.set_xlabel('Size'); ax.set_ylabel('Sphericity')
ax.set_title('Sphericity vs size', fontweight='bold')
ax.legend(fontsize=8); ax.grid(alpha=0.3)

# Inter-clone
ax = fig.add_subplot(gs[1, 0])
if len(inter_d) > 0:
    ax.hist(inter_d/1000, bins=30, color='red', alpha=0.7, edgecolor='black')
    ax.axvline(np.median(inter_d)/1000, color='blue', linestyle='--',
               label=f'med={np.median(inter_d)/1000:.2f}mm')
    ax.set_xlabel('Inter-clone distance (mm)'); ax.set_ylabel('Count')
    ax.set_title(f'{len(inter_d)} pairs', fontweight='bold')
    ax.legend(fontsize=8); ax.grid(alpha=0.3)

# 3D scatter
ax = fig.add_subplot(gs[1, 1], projection='3d')
nm = labels == -1
ax.scatter(centroids[nm,2], centroids[nm,1], centroids[nm,0], s=1, alpha=0.1, c='gray')
clone_cmap = plt.cm.tab20
for i, f in enumerate(clone_features[:20]):
    cm = labels == f['cluster_id']
    color = clone_cmap(i % 20)
    ax.scatter(centroids[cm,2], centroids[cm,1], centroids[cm,0], s=15, alpha=0.8, c=[color])
ax.set_xlabel('X', fontsize=8); ax.set_ylabel('Y', fontsize=8); ax.set_zlabel('Z', fontsize=8)
ax.set_title('Top 20 clones in 3D', fontweight='bold')

# Summary text
ax = fig.add_subplot(gs[1, 2])
ax.axis('off')
total_clone_cells = sum(f['n_cells'] for f in clone_features)
sm = f"""Z < {Z_MAX_UM}µm CLONE ANALYSIS

(Uniform detection range only)

Detection:
  Cells (Z<{Z_MAX_UM}µm):   {n_total:,}
  ({100*n_total/n_before:.1f}% of full {n_before})
  3D hull volume:    {hull_vol_mm3:.2f} mm³
  Density:           {n_total/hull_vol_mm3:.0f}/mm³

Spatial pattern:
  Mean NN:           {mean_nn:.1f} µm
  Clark-Evans R:     {clark_evans:.3f}

Clustering (eps={EPS_UM}µm):
  Total clusters:    {n_clusters}
  Largest:           {sizes_arr[0]} cells

Monte Carlo null:
  Clone thr:         ≥{clone_threshold} cells

★ CLONE CANDIDATES:
  N clones:          {n_clones}
  Cells in clones:   {total_clone_cells:,}
                     ({100*total_clone_cells/n_total:.1f}%)
  Largest clone:     {clone_features[0]['n_cells'] if clone_features else 0}
  Top 10 sizes:      {[f['n_cells'] for f in clone_features[:10]]}

Inter-clone (median): {np.median(inter_d)/1000 if len(inter_d)>0 else 0:.2f} mm
"""
ax.text(0.02, 0.98, sm, transform=ax.transAxes, fontsize=9,
        verticalalignment='top', fontfamily='monospace')

plt.suptitle(f'BDC tumor (DoG, Z<{Z_MAX_UM}µm): clone analysis',
             fontsize=13, fontweight='bold')
plt.savefig(OUTPUT_PNG, dpi=120, bbox_inches='tight')
plt.savefig(OUTPUT_SVG, format='svg', bbox_inches='tight')
plt.close()

print()
print('=' * 70)
print('  Summary (Z < {} µm):'.format(Z_MAX_UM))
print('=' * 70)
print(f"""
Detection:
  Cells (Z<{Z_MAX_UM}µm):   {n_total:,} ({100*n_total/n_before:.1f}% of {n_before})
  Hull volume:       {hull_vol_mm3:.2f} mm³
  Density:           {n_total/hull_vol_mm3:.0f} cells/mm³

Clark-Evans R:       {clark_evans:.3f} (clustered)

Clustering:
  N clusters:        {n_clusters}
  Largest cluster:   {sizes_arr[0]} cells
  
Monte Carlo null 99%ile: {clone_threshold}
★ N clones (≥{clone_threshold}):     {n_clones}
   Largest clone:    {clone_features[0]['n_cells'] if clone_features else 0}
   Top 10:           {[f['n_cells'] for f in clone_features[:10]]}

Saved: {OUTPUT_PNG}
""")
