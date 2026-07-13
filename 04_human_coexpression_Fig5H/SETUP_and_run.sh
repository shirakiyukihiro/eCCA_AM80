# ============================================================
# Fig 5H — environment setup and run commands
# Study      : Distinct CAF lineages and an Am80 stromal-reprogramming
#              strategy in extrahepatic cholangiocarcinoma (Meflin/ISLR).
# Repository : eCCA_AM80  (https://github.com/shirakiyukihiro/eCCA_AM80)
# Conda env + example run commands (paths generalised).
# Paths made relative for public release: put inputs in ./data, outputs in ./output.
# ============================================================

# Fig 5H co-expression — environment setup (same style as the napari env)
#
# Builds a conda env "coexpr" inside the writable PyTorch/CUDA container, with
# Cellpose (PyTorch, GPU via --nv) for nuclei segmentation + the scientific
# stack. The analysis also runs WITHOUT Cellpose (pass --no-cellpose) using a
# classical watershed fallback, so you can test immediately if the install lags.

# ============================================================ 1) create the env
sudo singularity shell --writable <container.sif>
source ~/.bashrc
conda create -n coexpr -c conda-forge python=3.11 mamba -y
conda activate coexpr

# scientific stack (imagecodecs handles LZW-compressed Keyence TIFFs)
mamba install -c conda-forge numpy scipy scikit-image pandas matplotlib tifffile imagecodecs -y

# PyTorch matching the container's CUDA 11.8, then Cellpose.
# Pin Cellpose to a 3.x release so the classic 'nuclei' model + API are stable.
python -m pip install --no-cache-dir torch --index-url https://download.pytorch.org/whl/cu118
python -m pip install --no-cache-dir "cellpose>=3.0,<4.0"

# sanity: confirm GPU is visible (should print True)
python -c "import torch; print('torch', torch.__version__, 'cuda', torch.cuda.is_available())"
exit

# ============================================================ 2) run the analysis
singularity exec --nv \
    --bind /path/to/data:/data:rw,/tmp:/tmp:rw \
    --env DISPLAY=$DISPLAY \
    <container.sif> \
    bash
source ~/.bashrc
conda activate coexpr

# (a) PROTOTYPE on ONE field, eyeball the QC images, decide the threshold K:
python fig5h_coexpr.py \
  "data/HBDC_Meflin_Thy1/<field_folder>" \
  --out output/qc_one --k 3

# (b) BATCH the whole study (root with all patient field folders):
python fig5h_coexpr.py \
  "data/HBDC_Meflin_Thy1" \
  --out output/fig5h_all --k 3

# CPU-only / no Cellpose (works anywhere): add  --no-cellpose
# options: --k 3.0  --ring-um 2.0  --sweep 1.5,2,2.5,3,4,5

# ============================================================ 3) outputs
#   <out>/per_image_summary.csv     n_Meflin_only / n_Thy1_only / n_double per field (CP932)
#   <out>/per_cell.csv              every cell: coords, intensities, class
#   <out>/per_patient_summary.csv   per-patient mean of field proportions (patient = unit)
#   <out>/<field>__seg_overlay.png  segmentation + classified cells
#   <out>/<field>__scatter.png      threshold-FREE Meflin x Thy1 per-cell scatter
#   <out>/<field>__sensitivity.png  double%% vs threshold K (robustness)
#   <out>/per_patient_proportions.png  study-level bar (mean +/- s.e.m., n patients)
