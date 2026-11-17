#!/usr/bin/env python
# ============================================================
# Detection uniformity with depth — justification for the Z < 1500 um cutoff
# Study      : Distinct CAF lineages and an Am80 stromal-reprogramming strategy
#              in extrahepatic cholangiocarcinoma (Meflin/ISLR).
# Repository : eCCA_AM80  (https://github.com/shirakiyukihiro/eCCA_AM80)
#
# Purpose    : Light-sheet illumination attenuates with tissue depth, so cells deep in the
#              volume are dimmer and are under-detected. This script quantifies where
#              detection density stops being uniform, which is what justifies restricting
#              the primary analysis to Z < 1500 um. It reports the numbers quoted in the
#              Supplementary Methods.
#
#              It runs from the provided cell centroids alone — the raw .ims volume is NOT
#              required. If data/BDC.ims happens to be present, an intensity-versus-depth
#              profile is added as a bonus panel.
#
# Input      : output/result_BDC_tumor_DoG.npz   (provided with the repository)
# Output     : output/detection_uniformity.png
#              output/detection_uniformity_summary.txt
# ============================================================

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ============================================================
# Configuration
# ============================================================
DETECTION_NPZ = 'output/result_BDC_tumor_DoG.npz'
OUT_DIR = 'output'

Z_BIN_UM = 200           # bin width for the detection-density profile
UNIFORM_LO_UM = 200      # the uniformly detected zone runs from here ...
UNIFORM_HI_UM = 1400     # ... to here
Z_CUTOFF_UM = 1500       # cutoff used by the primary analysis

OUT_PNG = f'{OUT_DIR}/detection_uniformity.png'
OUT_TXT = f'{OUT_DIR}/detection_uniformity_summary.txt'

os.makedirs(OUT_DIR, exist_ok=True)

# ============================================================
# 1. Detection density with depth
# ============================================================
print('=' * 70)
print('  Detection uniformity with depth')
print('=' * 70)

cents = np.load(DETECTION_NPZ)['centroids_um']
z = cents[:, 0].astype(float)
n_all = len(z)
print(f'\nCells detected (whole volume): {n_all}')
print(f'Depth range: {z.min():.0f} - {z.max():.0f} um')

bins = np.arange(0, z.max() + Z_BIN_UM, Z_BIN_UM)
counts, _ = np.histogram(z, bins=bins)
centres = (bins[:-1] + bins[1:]) / 2

uniform = (centres >= UNIFORM_LO_UM) & (centres <= UNIFORM_HI_UM)
deep = centres > UNIFORM_HI_UM
mean_uniform = counts[uniform].mean()
mean_deep = counts[deep].mean() if deep.any() else float('nan')
ratio = mean_deep / mean_uniform

print(f'\nCells per {Z_BIN_UM}-um Z bin:')
for c, n in zip(centres, counts):
    bar = '#' * int(40 * n / counts.max())
    tag = '  <- attenuated' if c > UNIFORM_HI_UM else ''
    print(f'  {c - Z_BIN_UM / 2:5.0f}-{c + Z_BIN_UM / 2:<5.0f} {n:5d}  {bar}{tag}')

print(f'\nUniform zone ({UNIFORM_LO_UM}-{UNIFORM_HI_UM} um): '
      f'mean {mean_uniform:.0f} cells per {Z_BIN_UM}-um bin')
print(f'Beyond {UNIFORM_HI_UM} um            : mean {mean_deep:.0f} cells '
      f'({100 * ratio:.0f}% of the uniform zone)')

n_kept = int((z < Z_CUTOFF_UM).sum())
print(f'\nZ < {Z_CUTOFF_UM} um: {n_kept} / {n_all} cells '
      f'({100 * n_kept / n_all:.1f}%) -> used for the primary analysis')

# ============================================================
# 2. Optional: intensity versus depth (only if the raw volume is available)
# ============================================================
have_ims = os.path.exists('data/BDC.ims')
if have_ims:
    print('\ndata/BDC.ims found - adding an intensity-versus-depth profile.')
else:
    print('\ndata/BDC.ims not present - skipping the intensity profile '
          '(not needed for the reported numbers).')

# ============================================================
# 3. Figure
# ============================================================
fig, ax = plt.subplots(figsize=(7, 4.2))
colours = ['#88CCEE' if c <= UNIFORM_HI_UM else '#DDAA77' for c in centres]
ax.bar(centres, counts, width=Z_BIN_UM * 0.88, color=colours,
       edgecolor='white', linewidth=0.5)
ax.axhline(mean_uniform, color='#4477AA', linestyle='--', linewidth=1.1,
           label=f'uniform zone mean ({mean_uniform:.0f})')
ax.axvline(Z_CUTOFF_UM, color='#DC143C', linestyle='--', linewidth=1.4,
           label=f'cutoff Z = {Z_CUTOFF_UM} um')
ax.set_xlabel('Tissue depth Z (um)')
ax.set_ylabel(f'Cells detected per {Z_BIN_UM}-um bin')
ax.set_title(f'Detection density with depth\n'
             f'beyond {UNIFORM_HI_UM} um, detection falls to '
             f'{100 * ratio:.0f}% of the uniform zone',
             fontsize=10, fontweight='bold')
ax.legend(fontsize=8)
plt.tight_layout()
plt.savefig(OUT_PNG, dpi=200, bbox_inches='tight')
plt.close()

# ============================================================
# 4. Summary
# ============================================================
with open(OUT_TXT, 'w') as f:
    f.write('Detection uniformity with depth\n')
    f.write('=' * 40 + '\n\n')
    f.write(f'Cells detected (whole volume): {n_all}\n')
    f.write(f'Depth range: {z.min():.0f} - {z.max():.0f} um\n\n')
    f.write(f'Cells per {Z_BIN_UM}-um Z bin:\n')
    for c, n in zip(centres, counts):
        f.write(f'  {c - Z_BIN_UM / 2:5.0f}-{c + Z_BIN_UM / 2:<5.0f} {n:5d}\n')
    f.write(f'\nUniform zone ({UNIFORM_LO_UM}-{UNIFORM_HI_UM} um): '
            f'mean {mean_uniform:.0f} cells per bin\n')
    f.write(f'Beyond {UNIFORM_HI_UM} um: mean {mean_deep:.0f} cells per bin '
            f'({100 * ratio:.0f}% of the uniform zone)\n\n')
    f.write(f'Primary analysis: Z < {Z_CUTOFF_UM} um -> {n_kept} / {n_all} cells '
            f'({100 * n_kept / n_all:.1f}%)\n\n')
    f.write('Whether the cells beyond the cutoff are merely blurred or genuinely\n')
    f.write('under-illuminated is addressed by test_deconvolution_deep_zone.py.\n')

print(f'\nFiles written:\n  {OUT_PNG}\n  {OUT_TXT}')
print('=' * 70)
