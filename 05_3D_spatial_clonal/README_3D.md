# 3D spatial analysis of Meflin-lineage cells in extrahepatic cholangiocarcinoma

Light-sheet analysis pipeline for:

> Ito K, Shiraki Y, *et al.* Identification of distinct cancer-associated fibroblast
> lineages reveals a therapeutic strategy for extrahepatic cholangiocarcinoma.

This folder covers two things:

1. **Three-dimensional renderings** of CUBIC-cleared, lineage-traced tissue
   (Figures 2D, 2H, 4F and Videos S1–S3) — produced interactively; see
   [`napari_visualisation_workflow.md`](napari_visualisation_workflow.md).
2. **Quantitative spatial analysis** of Meflin-lineage (tdTomato⁺) cells in the
   orthotopic tumour — cell detection, clustering statistics, and the
   **Figure 4G** density map.

---

## Files in this folder

| File | Role |
|---|---|
| `detect_cells_DoG.py` | **(1)** Difference-of-Gaussians cell detection in `.ims` volumes, inside region masks |
| `analyze_tumor_clones_DoG_Zcrop.py` | **(2)** Clark–Evans index, DBSCAN clustering, Monte-Carlo CSR null (Z < 1500 µm) |
| `figure_4G_kde_map.py` | **(3)** Renders the published **Fig. 4G** panel |
| `_dbscan_fallback.py` | DBSCAN implementation used if scikit-learn is unavailable |
| `validate_detection_uniformity.py` | Validation: signal attenuation with depth; justifies the Z < 1500 µm cutoff |
| `test_deconvolution_middle_zone.py` | Validation: Richardson-Lucy deconvolution test inside the uniformly detected zone |
| `test_deconvolution_deep_zone.py` | Validation: same test beyond the Z < 1500 µm cutoff (light-attenuated region) |
| `napari_visualisation_workflow.md` | How the 3D renderings and videos were produced |
| `environment.yml`, `requirements.txt` | Python environment |

---

## Reported numbers (BDC tumour, DoG detection, Z < 1500 µm)

These are the values quoted in the Results text.

| Metric | Value |
|---|---|
| Cells detected | **4,251** |
| Analysed volume (convex hull) | **14.01 mm³** |
| Cell density | **303 cells mm⁻³** |
| Mean nearest-neighbour distance | **53.5 µm** |
| Expected under CSR | **82.4 µm** |
| **Clark–Evans R** | **0.649** (R < 1 → clustered) |
| DBSCAN clusters (ε = 50 µm) | **307** |
| Clusters above the Monte-Carlo null (≥5 cells) | **127** |
| Largest cluster | **75 cells** |

> **Note on the Clark–Evans index.** In three dimensions the expected
> nearest-neighbour distance under complete spatial randomness is
> `E[d] = Γ(4/3) · ((4/3)·π·ρ)^(−1/3) = 0.55396 · ρ^(−1/3)`.
> An earlier version of `analyze_tumor_clones_DoG_Zcrop.py` used the
> **two-dimensional** coefficient (0.5), which gave R = 0.719. The corrected
> three-dimensional calculation gives **R = 0.649**, i.e. the cells are somewhat
> *more* clustered than previously reported. The conclusion is unchanged.

The observed mean nearest-neighbour distance (53.5 µm) is close to the DBSCAN
neighbourhood radius (ε = 50 µm), which is why that value was used.

> **Terminology.** Some file names, variable names and output keys in this code use the word
> *clone* (e.g. `analyze_tumor_clones_DoG_Zcrop.py`, `clone_threshold`), which reflects the
> original working hypothesis. The manuscript deliberately refers to these groups as **spatial
> clusters**, not clones: spatial proximity alone does not establish clonality, which would
> require lineage barcoding. Read every occurrence of *clone* in the code as *spatial cluster*.

---

## Environment

```bash
conda env create -f environment.yml
conda activate napari
```

Or, in a Singularity container with bind mounts:

```bash
singularity exec --nv --bind /your/data:/data your_image.sif bash
conda activate napari
```

---

## Data files

Raw light-sheet `.ims` volumes are **not distributed** (large; available from the
corresponding author on request). The region masks needed to reproduce the analysis
are included in `data/mask_polygons/`:

| File | Region |
|---|---|
| `mask_polygons_BDC_level3.json` | BDC whole tissue |
| `mask_polygons_tumor_BDC.json` | BDC tumour (used for Fig. 4G) |
| `mask_polygons_normal_BDC.json` | BDC peritumoural normal |
| `mask_polygons_normal_CBD.json` | CBD distant normal |
| `mask_polygons_Thy1.json` | Thy1-CreERT2 (SLICK-H) comparison |

Polygon format: per-Z-slice 2D contours, vertices as `[z, y, x]` in µm.

---

## Running the quantitative analysis

All paths are relative: inputs in `data/`, results written to `output/`.

### 1. Cell detection

```bash
python detect_cells_DoG.py
```

Edit the `SAMPLES` dictionary at the top of the script to point at your `.ims` files.

Outputs, per sample:
- `output/result_<sample>_DoG.npz` — cell centroids (µm)
- `output/result_<sample>_DoG_meta.json` — detection parameters used

### 2. Spatial statistics and clustering

```bash
python analyze_tumor_clones_DoG_Zcrop.py
```

Computes the Clark–Evans index, runs DBSCAN, and derives the cluster-size threshold
from a Monte-Carlo complete-spatial-randomness null (100 iterations, cells
redistributed within the same convex hull).

Outputs:
- `output/BDC_tumor_DoG_Zcrop_clones.npz` — centroids + cluster labels
- `output/BDC_tumor_DoG_Zcrop_cluster_features.csv` — per-cluster features
- `output/BDC_tumor_DoG_Zcrop_clone_analysis.svg` — overview figure (internal)

### 3. Figure 4G

```bash
python figure_4G_kde_map.py
```

Renders the published panel: every detected cell as a point (XY projection of the
three-dimensional coordinates), plus two kernel-density contour lines.

Outputs:
- `output/Figure_4G_KDE_map.svg` / `.png`

**What Fig. 4G shows.** The kernel-density estimate is computed on the XY projection
(`scipy.stats.gaussian_kde`, scalar bandwidth factor 0.03) on a 250 × 250 grid spanning
the data extent plus a 50-µm margin. Contours are drawn at the **92nd** and **99th**
percentiles of the grid density values — that is, they enclose the **densest 8%** (pink)
and the **densest 1%** (red) of the mapped area. The DBSCAN clusters are *not* drawn on
the panel; they are reported as numbers in the Results text.

### 4. Validation

```bash
python validate_detection_uniformity.py    # signal attenuation with depth
python test_deconvolution_middle_zone.py   # inside the uniformly detected zone
python test_deconvolution_deep_zone.py     # beyond the Z < 1500 um cutoff
```

`validate_detection_uniformity.py` profiles the mean in-mask intensity against tissue
depth and reports the **bottom/top intensity ratio** (deepest vs shallowest decile),
together with the number of cells detected per 100-µm Z bin. This is the basis for the
Z < 1500 µm cutoff.

The two deconvolution scripts each extract a 400 × 700 × 700 µm sub-volume centred on the
largest cluster in their zone, deconvolve it with an anisotropic Gaussian PSF, re-run DoG
detection with identical parameters, and write a before/after montage plus a summary. They
establish that the attenuation at depth is genuine signal loss rather than blur.

---

## Key parameters

| Parameter | Value |
|---|---|
| Imaris resolution level | 1 (highest) |
| Cell radius | 7.5 µm |
| DoG σ_small / σ_large | r/√2 / r·√2 |
| Minimum peak distance | 12 µm |
| DoG intensity threshold | 99.7th percentile (within mask) |
| Z range | 0 – 1500 µm |
| DBSCAN ε | 50 µm |
| DBSCAN min_samples | 3 |
| Cluster-size threshold | ≥5 cells (99th percentile of Monte-Carlo null) |
| Monte-Carlo iterations | 100 |
| Deconvolution PSF (validation) | Gaussian, σ_z = 7.83 µm, σ_xy = 2.35 µm |
| Richardson–Lucy iterations | 10 |
| KDE bandwidth factor | 0.03 (`gaussian_kde` scalar `bw_method`) |
| KDE grid | 250 × 250, data extent + 50 µm |
| KDE contour levels | 92nd, 99th percentile of grid density |

---

## Methodology notes

**Detection method.** Three approaches were evaluated: peak detection
(`min_distance` = 15 µm), finer peak detection (8 µm), and DoG-based 3D blob
detection. Peak detection under-detected closely packed cells; finer peak
detection introduced noise. DoG detection was selected on the basis of visual
validation and noise robustness.

**Z range.** Whole-tumour uniformity checks (`validate_detection_uniformity.py`) showed
light attenuation with depth (bottom/top intensity ratio ≈ 0.49). The primary analysis was
therefore restricted to Z < 1500 µm, where detection density is uniform.

**Deconvolution robustness.** Richardson–Lucy deconvolution (Gaussian point-spread
function, σ_xy = 2.35 µm, σ_z = 7.83 µm; 10 iterations) was tested on representative
sub-volumes with `test_deconvolution_middle_zone.py` and `test_deconvolution_deep_zone.py`.
It recovered ~5% additional cells in the middle zone (ratio 1.05×) and ~19% in the
deepest sub-volume (1.19×) — not enough to alter the primary results — so full-volume
deconvolution was not adopted. The modest gain at depth indicates that the attenuation
is genuine signal loss rather than blur, which is why the Z < 1500 µm cutoff is used.

**Single specimen.** The spatial analysis was performed on one representative tumour
and is reported descriptively; no between-animal inference is made.

---

## License

See `LICENSE`.

## Contact

Corresponding author (see manuscript).
