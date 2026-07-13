# 3D Spatial Analysis of Meflin⁺ Cells in Bile Duct Cancer

Quantitative analysis pipeline for:

> [Author et al., Title, Journal, Year, DOI]

This repository contains the workflow used to generate the figures in
the manuscript:
1. **napari 3D rendering screenshots** of Meflin-lineage-traced
   tdTomato⁺ cells in BDC tumor, peritumor normal, CBD distant normal,
   and Thy1Cre comparison samples.
2. **Quantitative spatial analysis figure** (Figure ?, panel B) — KDE
   density map of clonal clusters in BDC tumor.

## Quick Reference — Final Numbers

BDC tumor (DoG-based detection, Z < 1500 µm uniform zone):

| Metric                       | Value             |
|------------------------------|-------------------|
| Cells detected               | **4,251**         |
| Tumor volume analyzed        | **14.01 mm³**     |
| Cell density                 | **303 cells/mm³** |
| Clark–Evans R                | **0.719** (clustered) |
| DBSCAN clusters (ε = 50 µm)  | **307**           |
| Significant clones (≥5)      | **127**           |
| Largest clone                | **75 cells**      |
| Inter-clone median distance  | **~1.5 mm**       |

## Environment

```bash
conda env create -f environment.yml
conda activate napari
```

Or in a Singularity container with bind mounts, e.g.:
```bash
singularity exec --nv --bind /your/data:/data your_image.sif bash
conda activate napari
```

## Data Files

Raw `.ims` files are not included (large; available upon request).
Provided in `data/`:

| File                                 | Used by               |
|--------------------------------------|-----------------------|
| `mask_polygons_BDC_level3.json`      | BDC whole-tissue mask |
| `mask_polygons_tumor_BDC.json`       | BDC tumor (analysis)  |
| `mask_polygons_normal_BDC.json`      | BDC peritumor normal  |
| `mask_polygons_normal_CBD.json`      | CBD distant normal    |
| `mask_polygons_Thy1.json`            | Thy1Cre comparison    |

Polygon format: per-Z-slice 2D contours, vertices as `[z, y, x]` (µm).

---

## Section 1 — napari Visualization (image panels)

Used interactively to produce all 3D screenshots in the manuscript.

### Typical workflow

```bash
# Verify .ims voxel scale
python 01_napari_workflow/verify_ims_scale.py

# Open whole sample at level 3
python 01_napari_workflow/1_open_level3_full.py /path/to/sample.ims

# In napari console:
keep_inside = True
polygons_path = 'data/mask_polygons_tumor_BDC.json'
exec(open('01_napari_workflow/6_apply_saved_mask.py').read())

# Adjust contrast, take screenshot via File → Save Screenshot
# Save view state for reproducibility:
save_path = '/tmp/my_view.json'
exec(open('01_napari_workflow/8_save_view_state.py').read())
```

### Script reference

| Script                            | Purpose                          |
|-----------------------------------|----------------------------------|
| `verify_ims_scale.py`             | Sanity-check .ims voxel sizes    |
| `1_open_level3_full.py`           | Whole-sample view (level 3)      |
| `2_open_for_masking.py`           | Open for polygon annotation      |
| `3_save_mask_polygons.py`         | Save Shapes layer to JSON        |
| `4_load_mask_polygons.py`         | Load polygons into Shapes layer  |
| `5_apply_drawn_mask.py`           | Apply current Shapes layer mask  |
| `6_apply_saved_mask.py`           | Apply saved JSON mask            |
| `7_open_level0_crop.py`           | Open level-0 crop (high-res)     |
| `8_save_view_state.py`            | Save camera & contrast settings  |
| `9_reload_view_state.py`          | Restore camera & contrast        |
| `10_extract_oblique_slice.py`     | Tilted slice / MIP extraction    |
| `11_set_slab_clipping.py`         | 3D slab rendering                |
| `12_isolate_single_cell.py`       | Single-cell crop and rotation    |
| `get_view_crop_coords.py`         | Compute crop coordinates (µm)    |
| `inspect_auto_threshold.py`       | Threshold inspection             |
| `load_cluster_points.py`          | Overlay cluster centroids on the view |

---

## Section 2 — Quantitative Analysis

### 2.1 Cell detection (DoG)

```bash
python 02_quantitative_analysis/detect_cells_DoG_v3.py
```

Edit `SAMPLES` dictionary in script for your `.ims` paths.

Outputs (per sample):
- `result_<sample>_DoG.npz` — cell centroids in µm
- `result_<sample>_DoG_meta.json` — detection parameters used

### 2.2 Clone analysis (final, Z < 1500 µm)

```bash
python 02_quantitative_analysis/analyze_tumor_clones_DoG_Zcrop.py
```

Outputs:
- `BDC_tumor_DoG_Zcrop_clones.npz`           — cells + cluster labels
- `BDC_tumor_DoG_Zcrop_cluster_features.csv` — per-cluster features
- `BDC_tumor_DoG_Zcrop_clone_analysis.svg`   — overview figure

### 2.3 Main quantitative figure

```bash
python 04_figures/figure_kde_panelB_style.py
```

Outputs:
- `Figure_KDE_panelB_style.svg`
- `Figure_KDE_panelB_style.png`

---

## Section 3 — Validation (Methods documentation)

| Script                                  | What it documents             |
|-----------------------------------------|-------------------------------|
| `validate_detection_whole_tumor.py`     | Z-uniformity, intensity, top clones |
| `verify_DoG_peaks_montage_v2.py`        | Visual peak verification       |
| `visualize_new_largest_clone.py`        | Largest clone direct visualization |
| `test_deconvolution_subvolume.py`       | Deconvolution robustness (middle Z) |
| `test_deconvolution_deep_zone.py`       | Deconvolution robustness (deep Z) |

---

## Key Parameters

| Parameter                | Value             |
|--------------------------|-------------------|
| Imaris resolution level  | 1 (highest)       |
| Cell radius              | 7.5 µm            |
| DoG σ_small / σ_large    | r/√2  /  r·√2     |
| Min peak distance        | 12 µm             |
| DoG threshold            | 99.7th %ile (in-mask) |
| DBSCAN ε                 | 50 µm             |
| DBSCAN min_samples       | 3                 |
| Z range filter           | 0 – 1500 µm       |
| Clone threshold          | ≥ 5 cells (MC null 99%ile) |
| Monte Carlo iterations   | 100               |

## Methodology Notes

**Detection method.** We evaluated three approaches: peak detection
(`min_distance` = 15 µm), finer peak detection (8 µm), and DoG-based
3D blob detection. Peak detection under-detected closely-packed cells;
finer peak detection introduced noise; DoG-based detection was
selected as the final method based on visual validation and noise
robustness.

**Z range.** Whole-tumor uniformity validation revealed light
attenuation in deep tissue (intensity ratio bottom/top = 0.49).
Primary analysis was restricted to Z < 1500 µm where detection was
uniform.

**Deconvolution robustness.** Richardson–Lucy deconvolution was
tested on representative subvolumes. Middle-zone gain was 1.05× and
deep-zone gain was 1.19×, indicating that voxel sampling is the
resolution limit and deconvolution offers only marginal benefit. Full-
volume deconvolution was not adopted.

## License

[MIT recommended; see LICENSE]

## Citation

```
[Author et al., Title, Journal, Year, DOI]
```

## Contact

[Author email]
