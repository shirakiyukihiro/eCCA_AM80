# ============================================================
# PATCH for analyze_tumor_clones_DoG_Zcrop.py
#
# (1) BUG FIX — Clark-Evans expected nearest-neighbour distance
# (2) ADDITION — Monte-Carlo p-value for R (edge-corrected)
# ============================================================

# ------------------------------------------------------------
# (1) BUG FIX
# ------------------------------------------------------------
# The current code (line ~116) is:
#
#     expected_nn = 0.5 / (n_total / hull_vol_um3)**(1/3)
#
# The coefficient 0.5 is the *2-D* Clark-Evans constant (E[d] = 0.5 * rho**-0.5).
# In 3-D the expected nearest-neighbour distance under CSR is
#
#     E[d] = Gamma(4/3) / ((4/3) * pi * rho)**(1/3) = 0.55396 * rho**(-1/3)
#
# Using 0.5 makes the expected distance ~9.7% too small, which makes R ~10.8%
# too LARGE (i.e. it *under*-states the clustering). Replace with:

from math import gamma, pi

C_3D = gamma(4 / 3) / ((4 / 3) * pi) ** (1 / 3)   # = 0.553961...

rho = n_total / hull_vol_um3                       # cells per um^3
expected_nn = C_3D * rho ** (-1 / 3)
clark_evans = mean_nn / expected_nn
print(f'\nMean NN: {mean_nn:.1f}um, Expected CSR: {expected_nn:.1f}um, '
      f'CE R: {clark_evans:.3f}')

# Effect on the reported numbers:
#   n = 4,251 cells, hull volume = 14.01 mm^3 (rho = 303.4 cells/mm^3)
#   mean NN               = 53.5 um   (unchanged)
#   expected NN (old 0.5) = 74.4 um  -> R = 0.719   <- previously reported
#   expected NN (correct) = 82.4 um  -> R = 0.649   <- corrected value
#
# Note: R = 0.649 also means the observed mean NN (53.5 um) is close to the
# DBSCAN eps of 50 um, which is a defensible, data-driven justification for eps.


# ------------------------------------------------------------
# (2) ADDITION — Monte-Carlo p-value for R
# ------------------------------------------------------------
# The analytic expectation above assumes an unbounded domain. In a bounded convex
# hull, cells near the boundary have systematically longer nearest-neighbour
# distances (edge effect), which biases R upward. The Monte-Carlo loop already
# generates CSR realisations *inside the same hull*, so the null mean NN it
# produces is edge-corrected for free. Collect it and derive an empirical p-value.
#
# Inside the existing `for sim in range(N_MONTECARLO):` loop, right after
# `rand_pts = np.array(sampled[:n_total])`, add:

    rtree = cKDTree(rand_pts)
    rdd, _ = rtree.query(rand_pts, k=2)
    null_mean_nn.append(rdd[:, 1].mean())

# Before the loop, initialise:
#
#     null_mean_nn = []
#
# After the loop, add:

null_mean_nn = np.array(null_mean_nn)
R_mc = mean_nn / null_mean_nn.mean()               # edge-corrected R
n_more_extreme = int((null_mean_nn <= mean_nn).sum())
p_mc = (n_more_extreme + 1) / (N_MONTECARLO + 1)   # one-sided, add-one corrected

print(f'\nMonte-Carlo CSR null ({N_MONTECARLO} iterations):')
print(f'  null mean NN      : {null_mean_nn.mean():.1f} +/- '
      f'{null_mean_nn.std():.1f} um')
print(f'  observed mean NN  : {mean_nn:.1f} um')
print(f'  R (edge-corrected): {R_mc:.3f}')
print(f'  p (one-sided)     : {p_mc:.3f}'
      + ('  -> report as p < 0.01' if p_mc <= 0.01 else ''))

# With 100 iterations the smallest attainable p is 1/101 = 0.0099, so if the
# observed pattern is more clustered than every simulation you may report
# "p < 0.01 (100 Monte-Carlo iterations)". If you want a smaller p, raise
# N_MONTECARLO to 999 (-> p < 0.001).
