# ============================================================
# Fig 1D + 1E — Meflin/Thy1 fibroblast composition and co-expression
#               in biliary vs stromal area (mouse hepatoduodenal ligament)
#
# Study      : Distinct cancer-associated fibroblast lineages and an
#              Am80-based stromal-reprogramming strategy in extrahepatic
#              cholangiocarcinoma (Meflin/ISLR).
# Repository : eCCA_AM80  (https://github.com/shirakiyukihiro/eCCA_AM80)
#
# Purpose    : Read per-field Meflin (Islr) / Thy1 fluorescent-ISH counts in
#              the biliary and stromal compartments of the normal mouse
#              hepatoduodenal ligament, reduce each mouse to a single
#              proportion per compartment, and compare compartments with a
#              paired t-test (biological replicate = mouse, n = 5).
# Figure(s)  : Figure 1D (composition) and Figure 1E (co-expression).
# Input      : data/Fig1_data.csv.
# Output     : output/Fig1D.pdf, output/Fig1E.pdf.
# Style      : 00_plot_style.R (shared house style; box + beeswarm + p).
#
# Region mapping (verified against the figure):
#   n = Biliary (within ~75 um of bile-duct epithelium, "near")
#   d = Stromal (distal);  t = whole field -> excluded
#
# STATISTICS NOTE
#   The biological replicate is the MOUSE (n = 5), and each mouse provides
#   one biliary and one stromal value (paired). We therefore reduce each
#   mouse to a single proportion per region and compare regions with a
#   PAIRED t-test (n = 5 mice). We deliberately do NOT use a binomial GLMM
#   on pooled cell counts: that treats each of the hundreds of cells per
#   field as an independent trial, which is overdispersed / pseudoreplicated
#   and yields anti-conservative (falsely significant) p-values.
#
# Paths have been made relative for public release: place input files
# in ./data and write outputs to ./output (create these folders first).
# ============================================================

library(dplyr)
library(tidyr)
library(patchwork)
source("00_plot_style.R")

d0 <- read.csv("data/Fig1_data.csv", stringsAsFactors = FALSE)
d0 <- d0[d0$Region %in% c("d", "n"), ]
d0$Region <- factor(d0$Region, levels = c("n", "d"))   # Biliary, then Stromal

reg_lab     <- c(n = "Biliary", d = "Stromal")
cols_region <- c(n = "#18499e", d = "#e7211a")

## ---- per-mouse value (%) for the dots -----------------------------
## per-field proportion, then averaged within each mouse (one dot/mouse),
## consistent with the per-mouse averaging used for the IHC figures.
per_mouse <- function(num, den) {
  d0 %>%
    filter(!is.na(.data[[num]]), !is.na(.data[[den]]), .data[[den]] > 0) %>%
    mutate(p = .data[[num]] / .data[[den]] * 100) %>%
    group_by(MouseID, Region) %>%
    summarise(prop = mean(p), .groups = "drop")
}

## ---- paired t-test on the per-mouse proportions (n = 5 mice) -------
test_paired <- function(pm) {
  w <- pivot_wider(pm, names_from = Region, values_from = prop)
  w <- w[complete.cases(w[, c("n", "d")]), ]
  t.test(w$n, w$d, paired = TRUE)$p.value
}

panel <- function(num, den, ylab) {
  pm <- per_mouse(num, den)
  plot_box(pm, "Region", "prop", "Region",
           ylab = ylab, pval = test_paired(pm),
           cols = cols_region, xlabels = reg_lab,
           ylim = c(0, 100))                 # y-axis fixed 0-100 %
}

## ---- Fig 1D : composition (of total cells) ------------------------
pD1 <- panel("Thy1_Islr_Denom", "n_total",
             expression("Meflin"^"+" * " cells (%)"))
pD2 <- panel("Islr_Thy1_Denom", "n_total",
             expression("Thy1"^"+" * " cells (%)"))
fig1D <- pD1 + pD2 + plot_layout(ncol = 2)

## ---- Fig 1E : co-expression --------------------------------------
pE1 <- panel("Thy1_Islr_Num", "Thy1_Islr_Denom",
             expression("Thy1"^"+" * " in Meflin"^"+" * " cells (%)"))
pE2 <- panel("Islr_Thy1_Num", "Islr_Thy1_Denom",
             expression("Meflin"^"+" * " in Thy1"^"+" * " cells (%)"))
fig1E <- pE1 + pE2 + plot_layout(ncol = 2)

## ---- diagnostics: print per-mouse values + p for every panel ------
## (so you can see exactly what is computed; expected: the two
##  composition/co-expression "NS" panels are well above 0.05)
report <- function(num, den, label) {
  pm <- per_mouse(num, den)
  w  <- pivot_wider(pm, names_from = Region, values_from = prop)
  np <- sum(complete.cases(w[, c("n", "d")]))
  cat("\n----", label, "----\n")
  print(as.data.frame(w))
  cat(sprintf("paired t-test p = %.4g  (n = %d mice)\n", test_paired(pm), np))
}
cat("\n================ Fig 1 p-values (paired t-test, per-mouse) ================\n")
report("Thy1_Islr_Denom", "n_total",          "Fig 1D  Meflin+ cells (%)")
report("Islr_Thy1_Denom", "n_total",          "Fig 1D  Thy1+ cells (%)")
report("Thy1_Islr_Num",   "Thy1_Islr_Denom",  "Fig 1E  Thy1+ in Meflin+ (%)")
report("Islr_Thy1_Num",   "Islr_Thy1_Denom",  "Fig 1E  Meflin+ in Thy1+ (%)")
cat("==========================================================================\n")

## ---- export (sizes matched to the .ai panels; adjust if needed) ---
ggsave("output/Fig1D.pdf", fig1D, width = 84, height = 48, units = "mm", device = cairo_pdf)
ggsave("output/Fig1E.pdf", fig1E, width = 84, height = 48, units = "mm", device = cairo_pdf)

fig1D
fig1E
