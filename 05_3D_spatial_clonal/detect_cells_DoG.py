#!/usr/bin/env python
# ============================================================
# 3D cell detection by Difference-of-Gaussians (DoG) — FINAL detection method
# Study      : Distinct CAF lineages and an Am80 stromal-reprogramming strategy
#              in extrahepatic cholangiocarcinoma (Meflin/ISLR).
# Repository : eCCA_AM80  (https://github.com/shirakiyukihiro/eCCA_AM80)
# Detects Meflin-lineage (tdTomato+) cells in light-sheet .ims volumes within region masks. Outputs cell centroids (.npz) + parameters (.json).
# Paths made relative for public release: raw .ims files in ./data (not distributed;
# available on request), region masks in ./data/mask_polygons, results in ./output.
# ============================================================

"""
DoG (Difference of Gaussians) based cell detection - v3

v3 fixes:
  - JSON parser: polygons_um キー、µm 座標を voxel 変換
  - Z slice interpolation: 隣接 annotated Z 間を fill
  - IMSReader を with なしで使用 (v2 fix 継続)
"""
import os, sys, json
import numpy as np
from collections import defaultdict
from scipy.ndimage import gaussian_filter
from skimage.feature import peak_local_max
from matplotlib.path import Path as MplPath

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

SAMPLES = {
    'BDC_tumor': {
        'ims':           'data/BDC.ims',
        'channel':       1,
        'polygons':      'data/mask_polygons/mask_polygons_tumor_BDC.json',
        'output_npz':    'output/result_BDC_tumor_DoG.npz',
        'output_meta':   'output/result_BDC_tumor_DoG_meta.json',
    },
    'BDC_normal': {
        'ims':           'data/BDC.ims',
        'channel':       1,
        'polygons':      'data/mask_polygons/mask_polygons_normal_BDC.json',
        'output_npz':    'output/result_BDC_normal_DoG.npz',
        'output_meta':   'output/result_BDC_normal_DoG_meta.json',
    },
    'CBD_normal': {
        'ims':           'data/CBD.ims',
        'channel':       1,
        'polygons':      'data/mask_polygons/mask_polygons_normal_CBD.json',
        'output_npz':    'output/result_CBD_normal_DoG.npz',
        'output_meta':   'output/result_CBD_normal_DoG_meta.json',
    },
}

CELL_RADIUS_UM       = 7.5
MIN_PEAK_DISTANCE_UM = 12
DOG_RELATIVE_THR_PCT = 99.7
DETECT_LEVEL = 1


def open_ims(ims_path, level):
    from imaris_ims_file_reader.ims import ims as IMSReader
    return IMSReader(ims_path, ResolutionLevelLock=level)


def close_ims(imd):
    if hasattr(imd, 'close'):
        try: imd.close()
        except Exception: pass


def build_region_mask_v3(polygon_json, shape_zyx, voxel_um):
    """
    polygons_um (µm 座標) を読み込んで、各 Z slice にマスクを描画
    隣接 annotated Z 間は union で fill (swept volume)
    """
    with open(polygon_json) as f:
        data = json.load(f)
    
    polygons_um = data.get('polygons_um', [])
    n_total = data.get('n_polygons', len(polygons_um))
    print(f'  Total polygons: {n_total}')
    
    Z, Y, X = shape_zyx
    mask3d = np.zeros((Z, Y, X), dtype=bool)
    
    # Group polygons by their Z slice (Z は polygon 内で一定)
    polys_by_z = defaultdict(list)
    for poly in polygons_um:
        poly_arr = np.asarray(poly, dtype=float)
        if poly_arr.ndim != 2 or poly_arr.shape[1] != 3:
            continue
        # All Z 同じはず
        z_um = float(poly_arr[0, 0])
        z_idx = int(round(z_um / voxel_um[0]))
        
        # Convert (y_µm, x_µm) → voxel coords
        y_vox = poly_arr[:, 1] / voxel_um[1]
        x_vox = poly_arr[:, 2] / voxel_um[2]
        yx_vox = np.column_stack([y_vox, x_vox])
        polys_by_z[z_idx].append(yx_vox)
    
    z_keys = sorted(polys_by_z.keys())
    if z_keys:
        print(f'  Annotated Z slices: {len(z_keys)} unique')
        print(f'  Z slice range: {z_keys[0]} to {z_keys[-1]} (image Z = {Z})')
        # 各 Z での polygon 数
        z_counts = [(z, len(polys_by_z[z])) for z in z_keys]
        print(f'  Polygons per slice (top 5): {z_counts[:5]}')
        print(f'                       (last 5): {z_counts[-5:]}')
    
    # Pre-build coordinate grid
    yy, xx = np.meshgrid(np.arange(Y), np.arange(X), indexing='ij')
    pts_2d = np.column_stack([yy.ravel(), xx.ravel()])
    
    # Step 1: 各 annotated Z slice にマスク描画
    print(f'  Drawing polygons on annotated Z slices...')
    for z_idx, vert_list in polys_by_z.items():
        if z_idx < 0 or z_idx >= Z:
            continue
        slice_mask = np.zeros(Y * X, dtype=bool)
        for verts in vert_list:
            path = MplPath(verts)
            inside = path.contains_points(pts_2d)
            slice_mask |= inside
        mask3d[z_idx] = slice_mask.reshape(Y, X)
    
    annotated_voxels = mask3d.sum()
    print(f'  Voxels in annotated slices: {annotated_voxels:,}')
    
    # Step 2: 隣接 annotated Z 間を fill (swept volume = union of endpoints)
    print(f'  Filling Z gaps between annotated slices (swept union)...')
    annotated_z = sorted(polys_by_z.keys())
    if len(annotated_z) >= 2:
        for i in range(len(annotated_z) - 1):
            z1 = annotated_z[i]
            z2 = annotated_z[i + 1]
            if z2 - z1 <= 1:
                continue
            # Union of both endpoint masks
            union = mask3d[z1] | mask3d[z2]
            # Fill all intermediate slices
            for z in range(z1 + 1, z2):
                mask3d[z] = union
    
    # Step 3: 最初/最後の annotated Z の前後を少し extend (±2 slices)
    if annotated_z:
        first_z = annotated_z[0]
        last_z = annotated_z[-1]
        for z in range(max(0, first_z - 2), first_z):
            mask3d[z] |= mask3d[first_z]
        for z in range(last_z + 1, min(Z, last_z + 3)):
            mask3d[z] |= mask3d[last_z]
    
    total_voxels = mask3d.sum()
    vol_um3 = total_voxels * voxel_um[0] * voxel_um[1] * voxel_um[2]
    print(f'  Final mask voxels: {total_voxels:,}')
    print(f'  Mask volume: {vol_um3:,.0f} µm³ ({vol_um3/1e9:.3f} mm³)')
    
    return mask3d


def detect_cells_dog(ims_path, channel, level, polygon_json,
                       cell_radius_um, min_peak_dist_um, dog_thr_pct,
                       output_npz, output_meta, label='detect'):
    print(f'\n{"="*70}')
    print(f'  {label}: DoG-based detection (v3)')
    print(f'{"="*70}')
    
    # [1] Read .ims
    print(f'\n[1] Reading {ims_path} (level={level})...')
    imd = open_ims(ims_path, level)
    try:
        voxel_um = imd.resolution
        full_shape = imd.shape
        print(f'  Voxel: {voxel_um} µm')
        print(f'  Full shape: {full_shape}')
        
        if len(full_shape) == 5:
            img = np.asarray(imd[0, channel, :, :, :])
        elif len(full_shape) == 4:
            img = np.asarray(imd[channel, :, :, :])
        else:
            img = np.asarray(imd[:, :, :])
    finally:
        close_ims(imd)
    
    print(f'  Loaded: shape={img.shape}, dtype={img.dtype}, size={img.nbytes/1e9:.2f}GB')
    
    # [2] Build mask (v3 parser)
    print(f'\n[2] Building mask from {polygon_json}...')
    mask = build_region_mask_v3(polygon_json, img.shape, voxel_um)
    
    if mask.sum() == 0:
        print(f'\nERROR: Mask is still empty! Cannot proceed.')
        return None, None
    
    img_masked = img.astype(np.float32)
    img_masked[~mask] = 0
    del img
    
    # [3] DoG
    sigma_small_um = cell_radius_um / np.sqrt(2)
    sigma_large_um = cell_radius_um * np.sqrt(2)
    sigma_small_vox = tuple(sigma_small_um / float(v) for v in voxel_um)
    sigma_large_vox = tuple(sigma_large_um / float(v) for v in voxel_um)
    
    print(f'\n[3] DoG: σ_small={sigma_small_um:.1f}µm = {[f"{s:.2f}" for s in sigma_small_vox]} vox')
    print(f'        σ_large={sigma_large_um:.1f}µm = {[f"{s:.2f}" for s in sigma_large_vox]} vox')
    
    print(f'  Computing G(σ_small)...')
    g_small = gaussian_filter(img_masked, sigma_small_vox, mode='constant', cval=0)
    
    print(f'  Computing G(σ_large)...')
    g_large = gaussian_filter(img_masked, sigma_large_vox, mode='constant', cval=0)
    del img_masked
    
    print(f'  DoG = G_small - G_large')
    dog = g_small - g_large
    del g_small, g_large
    
    dog[~mask] = 0
    print(f'  DoG stats: min={dog.min():.2f}, max={dog.max():.2f}, mean={dog.mean():.4f}')
    
    # [4] Threshold
    dog_in_mask = dog[mask]
    if len(dog_in_mask) == 0:
        print(f'ERROR: No voxels in mask')
        return None, None
    
    dog_thr = float(np.percentile(dog_in_mask, dog_thr_pct))
    print(f'\n[4] DoG threshold ({dog_thr_pct}%ile in mask): {dog_thr:.2f}')
    n_above = (dog > dog_thr).sum()
    print(f'  Voxels above threshold: {n_above:,} ({100*n_above/dog.size:.4f}%)')
    
    # [5] Peaks
    min_dist_vox = max(1, int(min_peak_dist_um / min(voxel_um)))
    print(f'\n[5] peak_local_max (min_distance={min_dist_vox} vox = {min_peak_dist_um}µm)...')
    
    peaks_vox = peak_local_max(
        dog,
        min_distance=min_dist_vox,
        threshold_abs=dog_thr,
    )
    print(f'  Detected: {len(peaks_vox)} peaks')
    
    if len(peaks_vox) > 0:
        # Filter peaks within mask (peak_local_max should respect, but double check)
        in_mask = mask[peaks_vox[:, 0], peaks_vox[:, 1], peaks_vox[:, 2]]
        peaks_vox = peaks_vox[in_mask]
        print(f'  After mask filter: {len(peaks_vox)} peaks')
    
    del dog
    
    # [6] Save NPZ
    peaks_um = peaks_vox * np.array(voxel_um, dtype=np.float32)
    cluster_labels = -1 * np.ones(len(peaks_um), dtype=int)
    
    np.savez(
        output_npz,
        centroids_um=peaks_um.astype(np.float32),
        centroids_vox=peaks_vox.astype(np.int32),
        cluster_labels=cluster_labels,
    )
    print(f'\n[6] Saved NPZ: {output_npz}')
    
    meta = {
        'label': label, 'method': 'DoG v3',
        'ims_path': ims_path, 'level': level, 'channel': channel,
        'voxel_um': [float(v) for v in voxel_um],
        'cell_radius_um': cell_radius_um,
        'sigma_small_um': sigma_small_um,
        'sigma_large_um': sigma_large_um,
        'min_peak_distance_um': min_peak_dist_um,
        'dog_threshold_pct': dog_thr_pct,
        'dog_threshold_value': dog_thr,
        'n_cells': int(len(peaks_um)),
        'mask_voxels': int(mask.sum()),
        'mask_volume_um3': float(mask.sum() * voxel_um[0] * voxel_um[1] * voxel_um[2]),
    }
    with open(output_meta, 'w') as f:
        json.dump(meta, f, indent=2)
    print(f'  Saved meta: {output_meta}')
    
    print(f'\n========== Summary ({label}) ==========')
    print(f'  Mask volume: {meta["mask_volume_um3"]/1e9:.3f} mm³')
    print(f'  N cells: {len(peaks_um)}')
    print(f'  Density: {len(peaks_um) / (meta["mask_volume_um3"]/1e9):.0f} cells/mm³')
    
    return peaks_um, meta


if __name__ == '__main__':
    print('='*70)
    print('  DoG-based cell detection v3 (with proper polygon parser)')
    print(f'  r={CELL_RADIUS_UM}µm, min_dist={MIN_PEAK_DISTANCE_UM}µm, thr={DOG_RELATIVE_THR_PCT}%ile')
    print('='*70)
    
    for label, cfg in SAMPLES.items():
        try:
            detect_cells_dog(
                ims_path=cfg['ims'],
                channel=cfg['channel'],
                level=DETECT_LEVEL,
                polygon_json=cfg['polygons'],
                cell_radius_um=CELL_RADIUS_UM,
                min_peak_dist_um=MIN_PEAK_DISTANCE_UM,
                dog_thr_pct=DOG_RELATIVE_THR_PCT,
                output_npz=cfg['output_npz'],
                output_meta=cfg['output_meta'],
                label=label,
            )
        except Exception as e:
            import traceback
            print(f'\nERROR for {label}: {e}')
            traceback.print_exc()
    
    print('\n' + '='*70)
    print('  DONE. Next: python verify_DoG_peaks_montage_v2.py')
    print('='*70)
