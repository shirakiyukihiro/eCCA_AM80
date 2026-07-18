# eCCA_AM80 — analysis code

Analysis code for the study *"Identification of distinct cancer-associated fibroblast
lineages reveals a therapeutic strategy for extrahepatic cholangiocarcinoma"*
(Meflin/ISLR biliary fibroblasts; Am80/tamibarotene stromal reprogramming in
extrahepatic cholangiocarcinoma, eCCA).

## Repository layout

```
eCCA_AM80/
├── 00_plot_style.R              R — shared house style (box + beeswarm + p) for quantification figures
├── 01_scRNAseq/                 R — single-cell RNA-seq (mouse & human dCCA)
│   ├── 01_download_GSE163777.R                  Download mouse scRNA-seq (GEO GSE163777)
│   ├── 02_mouse_fibroblast_analysis.R           Mouse fibroblast Seurat + pseudotime (Figure 5A-E)
│   ├── 03_human_dCCA_annotation.R               Human dCCA cell-type annotation
│   ├── 04_human_dCCA_export_preCellBender_h5.R  Export .h5 for CellBender
│   ├── 04b_cellbender_remove_background.sh      CellBender ambient-RNA removal (per-sample params)
│   ├── 05_human_dCCA_postCellBender.R           Post-CellBender import & processing
│   └── 06_human_dCCA_figures.R                  Human dCCA UMAP / feature / violin figures (Figure 5F-H)
├── 02_image_quantification/     R — histology / ISH / IHC / morphometry quantification
│   ├── Fig1_DE.R                                Figure 1D/1E: Meflin/Thy1 composition & co-expression (paired t-test)
│   ├── clearing_morphometry.R                   Figure 2E: 3D-clearing morphometry, elongation/sphericity (Mann-Whitney)
│   ├── ISH_Islr_Thy1.R                          Figure 3B/D: ISH Islr/Thy1 (normal vs bile duct ligation; LMM)
│   ├── lineage_tracing_Cre_double.R             Figure 3H: Cre double-labelling co-expression (descriptive)
│   ├── IoTB_lineage_tracing.R                   Figure 4E: orthotopic (IoTB) lineage-tracing co-expression (descriptive)
│   ├── IHC_subcutaneous_WT_vs_KO.R              Figure 6E-H: subcutaneous IHC, WT vs Meflin-/- (LMM)
│   ├── IHC_subcutaneous_Am80_chemo.R            Figure 7E-J: subcutaneous IHC, Am80 + chemotherapy (LMM)
│   └── IHC_subcutaneous_Am80_PDL1.R             Figure 8C-L: subcutaneous IHC, Am80 + anti-PD-L1 (LMM)
├── 03_invivo_tumour/            R — in vivo tumour growth / body weight / endpoint volume
│   ├── Fig6_subcutaneous_WT_vs_KO.R             Figure 6B/C/D
│   ├── Fig7_Am80_chemotherapy.R                 Figure 7B/C/D
│   └── Fig8_S4_Am80_antiPDL1.R                  Figure 8B + Fig. S4A/B
├── 04_human_coexpression_Fig5H/ Python — Meflin/Thy1 co-expression in human dCCA
│   ├── fig5h_coexpr.py                          Cellpose segmentation + per-cell co-expression
│   ├── plot_fig5h.py                            Co-expression bar / pixel-colocalisation plots
│   ├── inspect_bcf_channels.py                  Keyence .bcf channel QC
│   └── SETUP_and_run.sh                         Environment setup + example commands
├── 05_3D_spatial_clonal/        Python — light-sheet 3D clonal spatial analysis
│   ├── detect_cells_DoG.py                      (1) DoG cell detection in .ims volumes
│   ├── analyze_tumor_clones_DoG_Zcrop.py        (2) Clark-Evans + DBSCAN + Monte-Carlo null
│   ├── _dbscan_fallback.py                      DBSCAN implementation (no scikit-learn needed)
│   ├── figure_4G_kde_map.py                     (3) renders the published Fig. 4G panel
│   ├── test_deconvolution_middle_zone.py        Validation: deconvolution test (uniform zone)
│   ├── test_deconvolution_deep_zone.py          Validation: deconvolution test (deep zone)
│   ├── README_3D.md                             Detailed 3D-analysis notes
│   └── environment.yml, requirements.txt        Python environment
├── data/
│   ├── Fig1_data.csv            Per-field Meflin/Thy1 ISH counts (Figure 1D/1E)
│   ├── *.xlsx                   Quantification workbooks (image / in vivo)
│   └── mask_polygons/           Region masks (JSON; per-Z polygons, world µm coordinates)
└── output/                      Figures and tables are written here
```

Folders 01–03 are **R**; folders 04–05 are **Python**. Each script begins with a header
block describing its purpose, the figure(s) it produces, and its input/output files.

## Running the code

All file paths are **relative**: put inputs in `data/` and write results to `output/`
(or edit the paths at the top of each script).

### scRNA-seq (01)
Scripts are numbered in execution order (01 → 06). Step **04b** is a shell script that runs
CellBender (GPU, Singularity) on the seven human dCCA samples; the per-sample parameters
reported in the manuscript are listed there.

### Image quantification / in vivo (02, 03)
Independent scripts; each reads one data file from `data/`. `Fig1_DE.R` reads
`data/Fig1_data.csv` and sources `00_plot_style.R` (shared box + beeswarm style).
`IHC_subcutaneous_WT_vs_KO.R`, `IHC_subcutaneous_Am80_chemo.R`, and
`IHC_subcutaneous_Am80_PDL1.R` share the same subcutaneous-IHC workbook and
analyse different comparison groups.

### 3D spatial / clonal analysis (05) — execution order
1. `detect_cells_DoG.py` — detects Meflin-lineage (tdTomato⁺) cells inside the region
   masks (`data/mask_polygons/`) and writes centroids to `output/`.
2. `analyze_tumor_clones_DoG_Zcrop.py` — computes the Clark-Evans index, clusters the
   centroids (DBSCAN, ε = 50 µm, min_samples = 3; Z < 1500 µm), and derives the
   cluster-size threshold from a Monte-Carlo complete-spatial-randomness null.
3. `figure_4G_kde_map.py` — renders the published **Fig. 4G** panel: all detected cells as
   points (XY projection), with 92nd- and 99th-percentile kernel-density contour lines
   enclosing the densest 8% and 1% of the mapped area. The DBSCAN clusters are reported as
   numbers in the Results text and are not drawn on the panel.
4. `test_deconvolution_middle_zone.py` / `test_deconvolution_deep_zone.py` — optional
   validation: Richardson-Lucy deconvolution (Gaussian PSF, σ_z = 7.83 µm, σ_xy = 2.35 µm,
   10 iterations) applied to a sub-volume, with DoG detection re-run before and after.
   These support the deconvolution-robustness statement in the Supplementary Methods.

Cell detection uses Difference-of-Gaussians (DoG) blob detection (cell radius 7.5 µm,
minimum peak distance 12 µm, 99.7th-percentile in-mask threshold), which was selected over
peak detection after visual validation.

Raw light-sheet `.ims` volumes are **not distributed** (large files; available on request).
The region masks required to reproduce the analysis are included in `data/mask_polygons/`.

## Software

R (Seurat, Slingshot, tradeSeq, clusterProfiler, Harmony, scDblFinder, lme4/lmerTest,
ggplot2/patchwork/ggbeeswarm) and Python (CellBender, Cellpose, napari, scikit-image,
scikit-learn, SciPy); see the manuscript's CTAT table and Supplementary Materials & methods
for full version numbers.

## Data availability

- Mouse bilio-vascular scRNA-seq: GEO **GSE163777**.
- Human dCCA scRNA-seq: obtained from the original authors on request (SRA **SUB11007007**).

## Notes

Some scripts contain working comments in Japanese; these do not affect execution.
