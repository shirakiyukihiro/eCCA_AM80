#!/usr/bin/env python
# ============================================================
# Deconvolution robustness test — DEEP zone (Z > 1500 um, light-attenuated region)
# Study      : Distinct CAF lineages and an Am80 stromal-reprogramming strategy
#              in extrahepatic cholangiocarcinoma (Meflin/ISLR).
# Repository : eCCA_AM80  (https://github.com/shirakiyukihiro/eCCA_AM80)
#
# Purpose    : The primary analysis is restricted to Z < 1500 um because light-sheet
#              signal attenuates with depth (bottom/top intensity ratio ~ 0.49). This
#              script asks whether the cells excluded by that cutoff are merely blurred
#              (recoverable by deconvolution) or genuinely under-illuminated. A
#              sub-volume centred on the largest deep cluster is deconvolved with the
#              same anisotropic Gaussian PSF used for the middle-zone test, and DoG
#              detection is re-run with identical parameters. Supports the statement in
#              the Supplementary Methods that deconvolution recovers ~19% additional
#              cells in the deepest sub-volume - not enough to alter the primary results,
#              so full-volume deconvolution was not adopted.
#
# NOTE       : Unlike the main pipeline, clusters here are computed on the FULL volume
#              (no Z cutoff), because the target lies beyond Z = 1500 um. DBSCAN is run
#              in-script with the same parameters as analyze_tumor_clones_DoG_Zcrop.py,
#              so only the detection output is needed as input.
#
# Input      : data/BDC.ims                            (not distributed; on request)
#              output/result_BDC_tumor_DoG.npz         (detect_cells_DoG.py)
# Output     : output/deconv_test_deep_montage.png
#              output/deconv_test_deep_summary.txt
# ============================================================

import os
import sys
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.ndimage import gaussian_filter
from skimage.restoration import richardson_lucy
from skimage.feature import peak_local_max
from collections import Counter

try:
    from sklearn.cluster import DBSCAN
    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False
    from _dbscan_fallback import dbscan_kdtree

# ============================================================
# Configuration
# ============================================================
IMS_PATH = 'data/BDC.ims'
DETECTION_NPZ = 'output/result_BDC_tumor_DoG.npz'
OUT_DIR = 'output'

CHANNEL = 1
LEVEL = 1

TARGET_Z_MIN_UM = 1500             # the depth beyond which the primary analysis stops
SUBVOL_HALF_UM = (200, 350, 350)   # same sub-volume size as the middle-zone test

# DBSCAN — identical to analyze_tumor_clones_DoG_Zcrop.py
EPS_UM = 50
MIN_SAMPLES = 3

NUM_ITER = 10                      # Richardson-Lucy iterations
PSF_SIGMA_VOX = (1.5, 1.0, 1.0)    # Gaussian PSF sigma in voxels, (Z, Y, X)
                                   # -> sigma_z = 7.83 um, sigma_xy = 2.35 um at level 1

# DoG detection parameters — identical to detect_cells_DoG.py
CELL_RADIUS_UM = 7.5
MIN_PEAK_DISTANCE_UM = 12
DOG_RELATIVE_THR_PCT = 99.7

MIDDLE_ZONE_RATIO = 1.05           # result of test_deconvolution_middle_zone.py

OUT_PNG = f'{OUT_DIR}/deconv_test_deep_montage.png'
OUT_TXT = f'{OUT_DIR}/deconv_test_deep_summary.txt'

os.makedirs(OUT_DIR, exist_ok=True)


# ============================================================
# Helpers
# ============================================================
def cluster_points(pts, eps, min_samples):
    if HAS_SKLEARN:
        return DBSCAN(eps=eps, min_samples=min_samples).fit(pts).labels_
    return dbscan_kdtree(pts, eps, min_samples)


def open_ims(path, level):
    from imaris_ims_file_reader.ims import ims as IMSReader
    return IMSReader(path, ResolutionLevelLock=level)


def close_ims(handle):
    if hasattr(handle, 'close'):
        try:
            handle.close()
        except Exception:
            pass


def dog_detect(image, voxel_um, cell_r_um, min_peak_um, thr_pct):
    """Difference-of-Gaussians blob detection — identical to the main pipeline."""
    sigma_s = cell_r_um / np.sqrt(2)
    sigma_l = cell_r_um * np.sqrt(2)
    sigma_s_vox = tuple(sigma_s / float(v) for v in voxel_um)
    sigma_l_vox = tuple(sigma_l / float(v) for v in voxel_um)

    img_f = image.astype(np.float32)
    g_s = gaussian_filter(img_f, sigma_s_vox, mode='constant', cval=0)
    g_l = gaussian_filter(img_f, sigma_l_vox, mode='constant', cval=0)
    dog = g_s - g_l
    del g_s, g_l, img_f

    thr = np.percentile(dog[dog != 0], thr_pct) if (dog != 0).sum() > 0 else 0
    min_dist_vox = max(1, int(min_peak_um / min(voxel_um)))
    peaks_vox = peak_local_max(dog, min_distance=min_dist_vox, threshold_abs=thr)
    return peaks_vox, thr


def gaussian_psf_3d(sigma_vox):
    """Normalised anisotropic 3D Gaussian PSF."""
    sig = np.asarray(sigma_vox)
    shape = tuple(max(3, int(6 * s + 1) | 1) for s in sig)
    coords = [np.arange(s) - s // 2 for s in shape]
    Z, Y, X = np.meshgrid(*coords, indexing='ij')
    g = np.exp(-(Z ** 2 / (2 * sig[0] ** 2)
                 + Y ** 2 / (2 * sig[1] ** 2)
                 + X ** 2 / (2 * sig[2] ** 2)))
    return g / g.sum()


# ============================================================
# 1. Cluster the FULL volume and pick the largest deep cluster
# ============================================================
print('=' * 70)
print(f'  Deconvolution robustness test — DEEP zone (Z > {TARGET_Z_MIN_UM} um)')
print('=' * 70)

print(f'\n[1] Loading detections from {DETECTION_NPZ} ...')
d = np.load(DETECTION_NPZ)
centroids = d['centroids_um']
print(f'  Cells detected in the full volume: {len(centroids)}')

print(f'\n[2] DBSCAN on the full volume (eps={EPS_UM} um, min_samples={MIN_SAMPLES}) ...')
labels = cluster_points(centroids, EPS_UM, MIN_SAMPLES)
counts = Counter(labels[labels != -1])
print(f'  Clusters: {len(counts)}')

deep = []
for cid, n in counts.items():
    cent = centroids[labels == cid].mean(axis=0)
    if cent[0] >= TARGET_Z_MIN_UM:
        deep.append((cid, n, cent))

if not deep:
    sys.exit(f'ERROR: no cluster with a centroid deeper than {TARGET_Z_MIN_UM} um.')

deep.sort(key=lambda t: -t[1])
print(f'\n  Deepest clusters (centroid Z > {TARGET_Z_MIN_UM} um):')
for cid, n, cent in deep[:5]:
    print(f'    id={cid:4d}  n={n:4d}  Z={cent[0]:.0f}  Y={cent[1]:.0f}  X={cent[2]:.0f}')

top_cid, top_n, center_um = deep[0]
print(f'\n  Selected: cluster {top_cid} ({top_n} cells) at Z = {center_um[0]:.0f} um')

# ============================================================
# 3. Read the sub-volume
# ============================================================
print('\n[3] Reading sub-volume from .ims ...')
imd = open_ims(IMS_PATH, LEVEL)
try:
    voxel = imd.resolution
    shape_full = imd.shape
    print(f'  Voxel size: {voxel} um')

    cz_v = int(center_um[0] / voxel[0])
    cy_v = int(center_um[1] / voxel[1])
    cx_v = int(center_um[2] / voxel[2])
    hz = max(1, int(SUBVOL_HALF_UM[0] / voxel[0]))
    hy = max(1, int(SUBVOL_HALF_UM[1] / voxel[1]))
    hx = max(1, int(SUBVOL_HALF_UM[2] / voxel[2]))

    if len(shape_full) == 5:
        _, _, Z, Y, X = shape_full
    elif len(shape_full) == 4:
        _, Z, Y, X = shape_full
    else:
        Z, Y, X = shape_full[-3:]

    zs, ze = max(0, cz_v - hz), min(Z, cz_v + hz + 1)
    ys, ye = max(0, cy_v - hy), min(Y, cy_v + hy + 1)
    xs, xe = max(0, cx_v - hx), min(X, cx_v + hx + 1)

    if len(shape_full) == 5:
        sub = np.asarray(imd[0, CHANNEL, zs:ze, ys:ye, xs:xe])
    elif len(shape_full) == 4:
        sub = np.asarray(imd[CHANNEL, zs:ze, ys:ye, xs:xe])
    else:
        sub = np.asarray(imd[zs:ze, ys:ye, xs:xe])
finally:
    close_ims(imd)

print(f'  Sub-volume: {sub.shape} ({sub.nbytes / 1e6:.1f} MB)')
print(f'  Z range:    {zs * voxel[0]:.0f} - {ze * voxel[0]:.0f} um')
print(f'  Intensity:  min {sub.min()}, median {np.median(sub):.0f}, '
      f'p95 {np.percentile(sub, 95):.0f}, max {sub.max()}')

# ============================================================
# 4-6. Detect -> deconvolve -> re-detect
# ============================================================
print('\n[4] DoG detection on the ORIGINAL sub-volume ...')
peaks_orig, thr_orig = dog_detect(
    sub, voxel, CELL_RADIUS_UM, MIN_PEAK_DISTANCE_UM, DOG_RELATIVE_THR_PCT)
print(f'  Cells detected: {len(peaks_orig)} (DoG threshold {thr_orig:.1f})')

print(f'\n[5] Richardson-Lucy deconvolution ({NUM_ITER} iterations) ...')
psf = gaussian_psf_3d(PSF_SIGMA_VOX)
print(f'  PSF in um: sigma_z = {PSF_SIGMA_VOX[0] * voxel[0]:.2f}, '
      f'sigma_y = {PSF_SIGMA_VOX[1] * voxel[1]:.2f}, '
      f'sigma_x = {PSF_SIGMA_VOX[2] * voxel[2]:.2f}')
t0 = time.time()
sub_norm = sub.astype(np.float32) / max(1, sub.max())
deconv = richardson_lucy(sub_norm, psf, num_iter=NUM_ITER, clip=False)
deconv_int = (deconv * sub.max()).astype(sub.dtype)
print(f'  Done in {time.time() - t0:.1f} s')

print('\n[6] DoG detection on the DECONVOLVED sub-volume ...')
peaks_deconv, thr_deconv = dog_detect(
    deconv_int, voxel, CELL_RADIUS_UM, MIN_PEAK_DISTANCE_UM, DOG_RELATIVE_THR_PCT)
print(f'  Cells detected: {len(peaks_deconv)} (DoG threshold {thr_deconv:.1f})')

ratio = len(peaks_deconv) / max(1, len(peaks_orig))
print('\n[7] Comparison')
print(f'  Original:    {len(peaks_orig):4d} cells')
print(f'  Deconvolved: {len(peaks_deconv):4d} cells')
print(f'  Ratio:       {ratio:.2f}x   (middle zone: {MIDDLE_ZONE_RATIO:.2f}x)')

# ============================================================
# 8. Montage
# ============================================================
print('\n[8] Writing montage ...')
N_PANELS = 8
z_indices = np.linspace(0, sub.shape[0] - 1, N_PANELS).astype(int)

fig = plt.figure(figsize=(N_PANELS * 2.5, 8))
gs = GridSpec(3, N_PANELS, figure=fig, hspace=0.18, wspace=0.05)

for row, (vol, peaks, name, colour) in enumerate([
        (sub, peaks_orig, 'ORIGINAL\n(DEEP)', 'cyan'),
        (deconv_int, peaks_deconv, 'DECONVOLVED\n(DEEP)', 'yellow')]):
    for i, z in enumerate(z_indices):
        ax = fig.add_subplot(gs[row, i])
        slc = vol[z]
        ax.imshow(slc, cmap='gray',
                  vmin=np.percentile(slc, 2), vmax=np.percentile(slc, 99.5))
        mask_z = np.abs(peaks[:, 0] - z) <= 2 if len(peaks) else None
        if mask_z is not None and mask_z.any():
            p = peaks[mask_z]
            ax.scatter(p[:, 2], p[:, 1], s=30, marker='*', c=colour,
                       edgecolor='red', linewidth=0.5)
        ax.set_title(f'Z={(z + zs) * voxel[0]:.0f}um\n'
                     f'({int(mask_z.sum()) if mask_z is not None else 0} cells +/-2)',
                     fontsize=8)
        ax.set_xticks([])
        ax.set_yticks([])
        if i == 0:
            ax.set_ylabel(name, fontsize=10, fontweight='bold')

for col, (vol, peaks, name, colour) in enumerate([
        (sub, peaks_orig, 'ORIGINAL (DEEP)', 'cyan'),
        (deconv_int, peaks_deconv, 'DECONVOLVED (DEEP)', 'yellow')]):
    ax = fig.add_subplot(gs[2, col * (N_PANELS // 2):(col + 1) * (N_PANELS // 2)])
    proj = vol.max(axis=0)
    ax.imshow(proj, cmap='gray',
              vmin=np.percentile(proj, 5), vmax=np.percentile(proj, 99.5))
    if len(peaks):
        ax.scatter(peaks[:, 2], peaks[:, 1], s=20, marker='*',
                   c=colour, edgecolor='red', linewidth=0.4)
    ax.set_title(f'{name}: max-Z projection ({len(peaks)} cells)',
                 fontsize=10, fontweight='bold')
    ax.set_xticks([])
    ax.set_yticks([])

fig.suptitle(f'Deconvolution test — DEEP zone (Z > {TARGET_Z_MIN_UM} um)\n'
             f'Cluster {top_cid} ({top_n} cells, Z = {center_um[0]:.0f} um) | '
             f'Gaussian PSF sigma_vox = {PSF_SIGMA_VOX}, '
             f'Richardson-Lucy {NUM_ITER} iter\n'
             f'Detection: {len(peaks_orig)} -> {len(peaks_deconv)} ({ratio:.2f}x); '
             f'middle zone {MIDDLE_ZONE_RATIO:.2f}x',
             fontsize=11, fontweight='bold')
plt.savefig(OUT_PNG, dpi=150, bbox_inches='tight')
plt.close()

# ============================================================
# 9. Summary
# ============================================================
with open(OUT_TXT, 'w') as f:
    f.write(f'Deconvolution robustness test - DEEP zone (Z > {TARGET_Z_MIN_UM} um)\n')
    f.write('=' * 55 + '\n\n')
    f.write(f'Target: largest deep cluster (id {top_cid}, {top_n} cells)\n')
    f.write(f'Centroid: Z={center_um[0]:.0f} Y={center_um[1]:.0f} '
            f'X={center_um[2]:.0f} um\n')
    f.write(f'Sub-volume: {sub.shape} voxels '
            f'(Z {zs * voxel[0]:.0f} - {ze * voxel[0]:.0f} um)\n\n')
    f.write('Intensity statistics (sub-volume):\n')
    f.write(f'  min {sub.min()}, median {np.median(sub):.0f}, '
            f'p95 {np.percentile(sub, 95):.0f}, max {sub.max()}\n\n')
    f.write(f'PSF (Gaussian): sigma = {PSF_SIGMA_VOX} voxels\n')
    f.write(f'                      = ({PSF_SIGMA_VOX[0] * voxel[0]:.2f}, '
            f'{PSF_SIGMA_VOX[1] * voxel[1]:.2f}, '
            f'{PSF_SIGMA_VOX[2] * voxel[2]:.2f}) um\n')
    f.write(f'Richardson-Lucy iterations: {NUM_ITER}\n\n')
    f.write('Detection (identical DoG parameters):\n')
    f.write(f'  Original:    {len(peaks_orig)} cells (threshold {thr_orig:.1f})\n')
    f.write(f'  Deconvolved: {len(peaks_deconv)} cells (threshold {thr_deconv:.1f})\n')
    f.write(f'  Ratio:       {ratio:.2f}x\n\n')
    f.write(f'Middle zone (test_deconvolution_middle_zone.py): '
            f'{MIDDLE_ZONE_RATIO:.2f}x\n\n')
    f.write('Interpretation: the gain at depth is modest, indicating that deep\n')
    f.write('attenuation reflects genuine signal loss rather than blur. Restricting\n')
    f.write(f'the primary analysis to Z < {TARGET_Z_MIN_UM} um therefore remains the\n')
    f.write('appropriate strategy, and full-volume deconvolution was not adopted.\n')

print(f'\nFiles written:\n  {OUT_PNG}\n  {OUT_TXT}')
print('=' * 70)
print(f'  DEEP zone:   {len(peaks_orig)} -> {len(peaks_deconv)} cells ({ratio:.2f}x)')
print(f'  MIDDLE zone: {MIDDLE_ZONE_RATIO:.2f}x')
print('=' * 70)
