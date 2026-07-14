# napari 3D visualisation workflow (light-sheet renderings)

This document records how the three-dimensional renderings of CUBIC-cleared, lineage-traced
tissue shown in the manuscript were generated. Rendering was performed **interactively in
napari**; the settings below are the ones used for the published images, so that the views can
be reproduced from the raw `.ims` volumes and the region masks in `data/mask_polygons/`.

Quantitative analysis of the same volumes (cell detection, clustering, Fig. 4G) is **not**
described here — see `detect_cells_DoG.py`, `analyze_tumor_clones_DoG_Zcrop.py`, and
`figure_kde_panelB_style.py`.

## Samples

| Volume | Lineage / tissue | Region masks |
|---|---|---|
| `Thy1Cre_stitch.ims` | Thy1-lineage (SLICK-H), extrahepatic bile duct | `mask_polygons_Thy1.json` |
| `CBD.ims` | Meflin-lineage, normal common bile duct | `mask_polygons_normal_CBD.json` |
| `BDC.ims` | Meflin-lineage, orthotopic bile duct tumour | `mask_polygons_BDC_level3.json` (whole tissue), `mask_polygons_tumor_BDC.json` (tumour), `mask_polygons_normal_BDC.json` (peritumoural normal) |

Channels: **tdTomato** (lineage-traced fibroblasts) and **CD31** (vessels).

## Environment

Python 3.11 with `napari`, `pyqt`, `napari-imaris-loader`, `aicsimageio`, `matplotlib`
(see `environment.yml`). Run with GPU/X11 access; `.ims` volumes are read through the
Imaris loader at a chosen resolution level (level 0 = full resolution).

## Steps

1. **Check voxel scale.** Confirm the µm voxel size reported by the `.ims` reader.

2. **Whole-sample overview (resolution level 3).** Open the volume, read µm coordinates from
   the status bar, and identify the region of interest.

3. **Region masks.** Exclusion/inclusion polygons were drawn on individual Z slices with the
   napari Shapes tool at a coarse level (4–5), interpolated between annotated slices, checked
   slice-by-slice, and saved as JSON (`data/mask_polygons/`, world µm coordinates). Masks were
   then re-applied at the level used for rendering. **The saved masks are provided, so this
   step does not need to be repeated.**

4. **Display settings.** Apply the saved mask, switch to 3D display, and set the contrast
   limits and rendering mode (below). Save the camera/contrast state so that the identical
   viewpoint can be restored on the high-resolution crop.

5. **High-resolution view (level 0 crop).** Re-open the region of interest at full resolution
   using the crop bounds below, re-apply the mask, restore the saved view state, adjust
   contrast, and capture the screenshot (scale bar 50 µm).

6. **Horizontal section** — zoomed-out view through the centre of the region (scale bar 200 µm).

7. **3D slab rendering.** Apply symmetric clipping planes about the view centre to render a
   slab of defined thickness, capture, then rotate the camera 90° about the view axis and
   capture again.

8. **Single-cell views.** Isolate individual tdTomato⁺ cells by restricting the volume to a
   sphere (radius 50–100 µm) centred on the cell, display the tdTomato channel only, zoom in,
   and capture (scale bar 10 µm).

## Settings used for the published images

| | `Thy1Cre_stitch.ims` | `CBD.ims` | `BDC.ims` |
|---|---|---|---|
| CD31 contrast limits | 1,500–15,000 | 1,500–15,000 | 500–2,000 |
| tdTomato contrast limits | 1,800–5,000 (level 3); 2,000–8,000 (level 0) | 5,000–40,000 | 2,726–10,601 |
| tdTomato rendering | `attenuated_mip` | `attenuated_mip` | `attenuated_mip` |
| Level-0 crop (z0 z1 y0 y1 x0 x1, µm) | 268 1742 1595 2638 1369 2985 | 500 2000 2704 3949 1090 2410 | 0 2200 3790 4390 50 3000 |
| Slab half-thickness | 50 µm (100 µm slab) | 50 µm (100 µm slab) | 100 µm (200 µm slab) |
| Single-cell isolation radius | 50–100 µm | 50–100 µm | — |

Scale bars: 200 µm (overview), 50 µm (high-resolution views), 10 µm (single cells).

## Note

The small interactive helper scripts used to drive napari (volume opening, mask application,
view-state save/restore, slab clipping, single-cell isolation) are available from the
corresponding author on request. All parameters needed to reproduce the published views are
listed above.
