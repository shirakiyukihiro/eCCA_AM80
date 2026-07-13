#!/usr/bin/env bash
# ============================================================
# Ambient-RNA removal with CellBender — human dCCA scRNA-seq
#
# Study      : Distinct cancer-associated fibroblast lineages and an Am80
#              stromal-reprogramming strategy in extrahepatic cholangiocarcinoma.
# Repository : eCCA_AM80  (https://github.com/shirakiyukihiro/eCCA_AM80)
#
# Purpose    : Removes ambient RNA from the seven human distal cholangiocarcinoma
#              (dCCA) samples, then repacks the corrected matrices so that they can
#              be read directly by Seurat.
# Pipeline   : 04_human_dCCA_export_preCellBender_h5.R  ->  THIS SCRIPT  ->
#              05_human_dCCA_postCellBender.R  ->  06_human_dCCA_figures.R
# Input      : data/human_dCCA_scRNAseq/matrix_1.h5 ... matrix_7.h5
# Output     : *_cellbender_filtered.h5  and  *_cellbender_filtered_seurat.h5
#
# Environment: CellBender v0.3.2, run in a Singularity container with GPU support
#              (docker://us.gcr.io/broad-dsde-methods/cellbender:latest).
#
# NOTE ON PARAMETERS
# ------------------
# Parameters were optimised PER SAMPLE by inspecting the CellBender QC report of a
# first, default run, following these rules:
#   --expected-cells            : cell number at the first steep drop of the UMI curve
#   --total-droplets-included   : cell number just past the plateau of the UMI curve
#   --learning-rate             : halved when the training/test ELBO curves diverged
#                                 or the test curve deteriorated
#   --epochs                    : increased (up to 300) when the ELBO had not converged
#   --fpr                       : raised above the (conservative) default of 0.01 when
#                                 ambient contamination was high
# The values below are the ones used for the results reported in the manuscript.
# Samples not listed with extra flags were run with CellBender defaults apart from --fpr.
# ============================================================
set -euo pipefail

DATA_DIR="data/human_dCCA_scRNAseq"

run_cellbender () {           # $1 = sample id, remaining args = sample-specific flags
    local S="$1"; shift
    local IN="${DATA_DIR}/matrix_${S}.h5"
    local OUT="${DATA_DIR}/matrix_${S}_cellbender.h5"

    echo "=== CellBender: sample ${S} ==="
    cellbender remove-background --cuda --input "${IN}" --output "${OUT}" "$@"

    # Repack the filtered matrix so that Seurat can read it directly
    ptrepack --complevel 5 \
        "${DATA_DIR}/matrix_${S}_cellbender_filtered.h5:/matrix" \
        "${DATA_DIR}/matrix_${S}_cellbender_filtered_seurat.h5:/matrix"
}

# ---- per-sample settings used for the manuscript --------------------------------
run_cellbender 1 --expected-cells 4000 --total-droplets-included 20000 \
                 --learning-rate 2.5e-5 --epochs 200 --fpr 0.05
run_cellbender 2 --fpr 0.05
run_cellbender 3 --expected-cells 3000 --total-droplets-included 40000 \
                 --learning-rate 1e-5   --epochs 200 --fpr 0.05
run_cellbender 4 --fpr 0.1
run_cellbender 5 --fpr 0.05
run_cellbender 6 --expected-cells 1500 --total-droplets-included 30000 \
                 --learning-rate 1e-5   --epochs 200 --fpr 0.05
run_cellbender 7 --fpr 0.1

echo "All samples processed."
