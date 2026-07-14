#!/usr/bin/env python
# ============================================================
# Detection uniformity with depth — justification for the Z < 1500 um cutoff
# Study      : Distinct CAF lineages and an Am80 stromal-reprogramming strategy
#              in extrahepatic cholangiocarcinoma (Meflin/ISLR).
# Repository : eCCA_AM80  (https://github.com/shirakiyukihiro/eCCA_AM80)
#
# Purpose    : Light-sheet illumination attenuates with tissue depth, so cells deep in
#              the volume are dimmer and are under-detected. This script quantifies that
#              attenuation inside the tumour mask and produces the two numbers used to
#              justify restricting the primary analysis to Z < 1500 um:
#
#                (a) bottom/top intensity ratio  — mean in-mask signal in the deepest
#                    Z decile divided by that in the shallowest Z decile;
#                (b) detected cells per Z bin    — the depth beyond which detection
#                    density falls away.
#
#              The Supplementary Methods quote a bottom/top intensity ratio of ~0.49.
#
# Input      : data/BDC.ims                                   (not distributed; on request)
#              data/mask_polygons/mask_polygons_tumor_BDC.json
#              output/result_BDC_tumor_DoG.npz                (detect_cells_DoG.py)
# Output     : output/detection_uniformity.png
#              output/detection_uniformity_summary.txt
# ============================================================

import os
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath

# ============================================================
# Configuration
# ============================================================
IMS_PATH = 'data/BDC.ims'
POLYGON_JSON = 'data/mask_polygons/mask_polygons_tumor_BDC.json'
DETECTION_NPZ = 'output/result_BDC_tumor_DoG.npz'
OUT_DIR = 'output'

CHANNEL = 1
LEVEL = 1

N_DECILES = 10        # the bottom/top ratio compares the deepest and shallowest decile
Z_BIN_UM = 100        # bin width for the detection-density profile
Z_CUTOFF_UM = 1500    # the cutoff used by the primary analysis

OUT_PNG = f'{OUT_DIR}/detection_uniformity.png'
OUT_TXT = f'{OUT_DIR}/detection_uniformity_summary.txt'

os.makedirs(OUT_DIR, exist_ok=True)


# ============================================================
# Helpers
# ============================================================
def open_ims(path, level):
    from imaris_ims_file_reader.ims import ims as IMSReader
    return IMSReader(path, ResolutionLevelLock=level)


def close_ims(handle):
    if hasattr(handle, 'close'):
        try:
            handle.close()
        except Exception:
            pass


def load_polygons_um(path):
    """Region mask JSON: per-Z 2D contours, vertices as [z, y, x] in micrometres."""
    with open(path) as f:
        data = json.load(f)
    polys = data['polygons_um'] if isinstance(data, dict) else data
    by_z = {}
    for poly in polys:
        arr = np.asarray(poly, dtype=float)
        z = float(np.mean(arr[:, 0]))
        by_z.setdefault(round(z), []).append(arr[:, 1:])   # (y, x) vertices
    return by_z


def mask_for_slice(polys_yx, shape_yx, voxel_yx):
    """Rasterise the (y, x) polygons of one Z slice onto a boolean mask."""
    yy, xx = np.mgrid[0:shape_yx[0], 0:shape_yx[1]]
    pts = np.column_stack([yy.ravel(), xx.ravel()])
    mask = np.zeros(pts.shape[0], dtype=bool)
    for poly_um in polys_yx:
        poly_vox = poly_um / np.asarray(voxel_yx)
        mask |= MplPath(poly_vox).contains_points(pts)
    return mask.reshape(shape_yx)


# ============================================================
# 1. Intensity versus depth, inside the tumour mask
# ============================================================
print('=' * 70)
print('  Detection uniformity with depth')
print('=' * 70)

print(f'\n[1] Loading region mask: {POLYGON_JSON}')
polys_by_z = load_polygons_um(POLYGON_JSON)
annotated_z = sorted(polys_by_z)
print(f'  Annotated Z slices: {len(annotated_z)} '
      f'({annotated_z[0]} - {annotated_z[-1]} um)')

print(f'\n[2] Reading {IMS_PATH} and profiling intensity with depth ...')
imd = open_ims(IMS_PATH, LEVEL)
try:
    voxel = imd.resolution                       # (z, y, x) um
    shape_full = imd.shape
    if len(shape_full) == 5:
        _, _, Z, Y, X = shape_full
    elif len(shape_full) == 4:
        _, Z, Y, X = shape_full
    else:
        Z, Y, X = shape_full[-3:]
    print(f'  Voxel {voxel} um, volume {Z} x {Y} x {X}')

    z_um_list, mean_int_list = [], []
    for z_key in annotated_z:
        z_vox = int(round(z_key / voxel[0]))
        if not (0 <= z_vox < Z):
            continue
        if len(shape_full) == 5:
            slc = np.asarray(imd[0, CHANNEL, z_vox, :, :])
        elif len(shape_full) == 4:
            slc = np.asarray(imd[CHANNEL, z_vox, :, :])
        else:
            slc = np.asarray(imd[z_vox, :, :])

        m = mask_for_slice(polys_by_z[z_key], slc.shape, (voxel[1], voxel[2]))
        if m.sum() == 0:
            continue
        z_um_list.append(z_key)
        mean_int_list.append(float(slc[m].mean()))
finally:
    close_ims(imd)

z_um = np.asarray(z_um_list, dtype=float)
mean_int = np.asarray(mean_int_list, dtype=float)
order = np.argsort(z_um)
z_um, mean_int = z_um[order], mean_int[order]
print(f'  Profiled {len(z_um)} slices inside the mask')

# ============================================================
# 2. Bottom/top intensity ratio
# ============================================================
k = max(1, len(z_um) // N_DECILES)
top_mean = mean_int[:k].mean()        # shallowest decile
bottom_mean = mean_int[-k:].mean()    # deepest decile
ratio = bottom_mean / top_mean

print('\n[3] Bottom/top intensity ratio')
print(f'  Shallowest decile (Z {z_um[0]:.0f}-{z_um[k - 1]:.0f} um): '
      f'mean intensity {top_mean:.1f}')
print(f'  Deepest decile    (Z {z_um[-k]:.0f}-{z_um[-1]:.0f} um): '
      f'mean intensity {bottom_mean:.1f}')
print(f'  bottom/top ratio = {ratio:.2f}')
print('  (the Supplementary Methods quote ~0.49)')

# ============================================================
# 3. Detected cells per Z bin
# ============================================================
print(f'\n[4] Detection density with depth ({DETECTION_NPZ})')
cents = np.load(DETECTION_NPZ)['centroids_um']
z_cells = cents[:, 0]
bins = np.arange(0, z_cells.max() + Z_BIN_UM, Z_BIN_UM)
counts, _ = np.histogram(z_cells, bins=bins)
centres = (bins[:-1] + bins[1:]) / 2

within = z_cells < Z_CUTOFF_UM
print(f'  Cells total          : {len(z_cells)}')
print(f'  Cells at Z < {Z_CUTOFF_UM} um : {int(within.sum())} '
      f'({100 * within.mean():.1f}%)')

# ============================================================
# 4. Figure
# ============================================================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))

ax1.plot(z_um, mean_int, color='#333333', linewidth=1.2)
ax1.axvline(Z_CUTOFF_UM, color='#DC143C', linestyle='--', linewidth=1.2,
            label=f'cutoff Z = {Z_CUTOFF_UM} um')
ax1.axhline(top_mean, color='#4477AA', linestyle=':', linewidth=1,
            label=f'shallowest decile ({top_mean:.0f})')
ax1.axhline(bottom_mean, color='#EE7733', linestyle=':', linewidth=1,
            label=f'deepest decile ({bottom_mean:.0f})')
ax1.set_xlabel('Tissue depth Z (um)')
ax1.set_ylabel('Mean in-mask intensity')
ax1.set_title(f'Signal attenuation with depth\nbottom/top ratio = {ratio:.2f}',
              fontsize=10, fontweight='bold')
ax1.legend(fontsize=8)

ax2.bar(centres, counts, width=Z_BIN_UM * 0.88,
        color='#88CCEE', edgecolor='white', linewidth=0.4)
ax2.axvline(Z_CUTOFF_UM, color='#DC143C', linestyle='--', linewidth=1.2,
            label=f'cutoff Z = {Z_CUTOFF_UM} um')
ax2.set_xlabel('Tissue depth Z (um)')
ax2.set_ylabel(f'Cells detected per {Z_BIN_UM}-um bin')
ax2.set_title('Detection density with depth', fontsize=10, fontweight='bold')
ax2.legend(fontsize=8)

plt.tight_layout()
plt.savefig(OUT_PNG, dpi=200, bbox_inches='tight')
plt.close()

# ============================================================
# 5. Summary
# ============================================================
with open(OUT_TXT, 'w') as f:
    f.write('Detection uniformity with depth\n')
    f.write('=' * 40 + '\n\n')
    f.write(f'Slices profiled inside the tumour mask: {len(z_um)}\n')
    f.write(f'Depth range: {z_um[0]:.0f} - {z_um[-1]:.0f} um\n\n')
    f.write('Intensity (mean signal inside the mask):\n')
    f.write(f'  shallowest decile (Z {z_um[0]:.0f}-{z_um[k - 1]:.0f} um): '
            f'{top_mean:.1f}\n')
    f.write(f'  deepest decile    (Z {z_um[-k]:.0f}-{z_um[-1]:.0f} um): '
            f'{bottom_mean:.1f}\n')
    f.write(f'  bottom/top ratio  : {ratio:.2f}\n\n')
    f.write('Detection:\n')
    f.write(f'  cells total            : {len(z_cells)}\n')
    f.write(f'  cells at Z < {Z_CUTOFF_UM} um   : {int(within.sum())} '
            f'({100 * within.mean():.1f}%)\n\n')
    f.write(f'The primary analysis is restricted to Z < {Z_CUTOFF_UM} um, where\n')
    f.write('detection density is uniform. Whether the cells beyond that depth are\n')
    f.write('merely blurred or genuinely under-illuminated is addressed by\n')
    f.write('test_deconvolution_deep_zone.py.\n')

print(f'\nFiles written:\n  {OUT_PNG}\n  {OUT_TXT}')
print('=' * 70)
print(f'  bottom/top intensity ratio = {ratio:.2f}')
print('=' * 70)
