# ============================================================
# Figure 8 / Supplementary S4 — Am80 + anti-PD-L1
#
# Study      : Distinct cancer-associated fibroblast lineages and an
#              Am80-based stromal-reprogramming strategy in extrahepatic
#              cholangiocarcinoma (Meflin/ISLR).
# Repository : eCCA_AM80  (https://github.com/shirakiyukihiro/eCCA_AM80)
#
# Purpose    : Endpoint tumour volume (one-way ANOVA + Tukey HSD) and tumour-growth / body-weight curves (mean ± SEM + pairwise LMM) for the Am80 + anti-PD-L1 experiment.
# Figure(s)  : Figure 8B and Supplementary Figure S4A/B.
# Input      : data/Am80_PDL1.xlsx.
# Output     : output/Fig8BS4AB.pdf.
#
# Paths have been made relative for public release: place input files
# in ./data and write outputs to ./output (create these folders, or edit
# the paths below).
# ============================================================

# ============================================================
# FIGURE 8  (Am80 + anti-PD-L1)  + Supplementary S4A/B
#   8B   Endpoint tumor volume -> one-way ANOVA + Tukey HSD (4 groups)
#   S4A  Tumor growth curves   -> mean ± SEM  + pairwise LMM (Am80+PD-L1 vs PD-L1)
#   S4B  Body weight curves     -> mean ± SEM  + pairwise LMM (Am80+PD-L1 vs PD-L1)
# Single genotype (WT). All error bars = SEM.
# ============================================================
# Data file: Am80_PDL1.xlsx
#   sheet 'Volume' : tumor volumes (col1 = group [first row of each group only]; wks in row 1)
#   sheet 'Sheet2' : body weights  (SAME layout; early cols wk3.5/4 hold TUMOR volumes -> filtered out)
#
# Treatment label -> Figure group:
#   DMSO+IgG  -> Control        Am80+IgG  -> Am80
#   DMSO+PDL1 -> PD-L1          Am80+PDL1 -> Am80 + PD-L1
# ------------------------------------------------------------

library(readxl)
library(tidyverse)
library(lme4)
library(lmerTest)
# resolve common namespace masks (Bioconductor etc.)
library(conflicted)
# dplyr / tidyr / lmerTest を優先
for (f in c("select","filter","mutate","rename","count","distinct","summarise",
            "summarize","group_by","arrange","transmute","slice")) conflict_prefer(f, "dplyr")
conflict_prefer("pivot_longer", "tidyr")
conflict_prefer("lmer", "lmerTest")
# base 関数を優先（Matrix / S4Vectors 等のマスクを回避）
for (f in c("unname","as.data.frame","colnames","rownames","setdiff","union",
            "intersect","which","table","cbind","rbind")) conflict_prefer(f, "base")

file_path <- "data/Am80_PDL1.xlsx"

grp_map8   <- c("DMSO+IgG"="Control","Am80+IgG"="Am80","DMSO+PDL1"="PD-L1","Am80+PDL1"="Am80 + PD-L1")
grp_lv8    <- c("Control","Am80","PD-L1","Am80 + PD-L1")
grp_cols8  <- c("Control"="grey55","Am80"="#377EB8","PD-L1"="#4DAF4A","Am80 + PD-L1"="#E41A1C")

## ---- read a sheet: fill-down group labels, pivot, return tidy long ----
read_fig8 <- function(path, sheet, value_name = "value", bw_filter = FALSE) {
  raw <- read_excel(path, sheet = sheet, col_names = FALSE)
  wk  <- suppressWarnings(as.numeric(as.character(unlist(raw[1, ]))))  # row1 = weeks (col1 = NA)
  body <- raw[-1, ]
  grp <- as.character(body[[1]])
  for (i in seq_along(grp)) if (i > 1 && (is.na(grp[i]) || grp[i] == "")) grp[i] <- grp[i-1]  # fill down
  
  meas <- body[, -1, drop = FALSE]
  meas <- dplyr::mutate(meas, dplyr::across(dplyr::everything(), as.character))  # unify type
  colnames(meas) <- as.character(wk[-1])
  meas$.gid   <- grp
  meas$.mouse <- paste0(sheet, "_", seq_len(nrow(meas)))
  
  long <- meas %>%
    pivot_longer(-c(.gid, .mouse), names_to = "week", values_to = "v") %>%
    dplyr::mutate(week  = suppressWarnings(as.numeric(week)),
                  value = suppressWarnings(as.numeric(v))) %>%
    dplyr::filter(!is.na(week), !is.na(value))
  if (bw_filter) long <- dplyr::filter(long, value < 40)   # drop tumor-volume contamination (wk3.5/4)
  long %>%
    dplyr::mutate(group = factor(unname(grp_map8[.gid]), levels = grp_lv8),
                  mouse = factor(.mouse)) %>%
    dplyr::transmute(mouse, group, week, !!value_name := value)
}

vol <- read_fig8(file_path, "Volume", "volume")
bw  <- read_fig8(file_path, "BW", "bw", bw_filter = TRUE)

## ---- SANITY CHECKS (expect 5 per group, 20 total) ----
cat("\n=== Fig8 mice per group (volume) — expect 5 each, 4 groups ===\n")
print(vol %>% dplyr::distinct(group, mouse) %>% dplyr::count(group))
cat("\n=== week range (volume) ===\n"); print(range(vol$week))
cat("\n=== body-weight value range (should be ~17-24 g, NO tumor volumes) ===\n"); print(range(bw$bw))
cat("=== body-weight week range ===\n"); print(range(bw$week))

# =================================================================
# 8B  ENDPOINT TUMOR VOLUME : one-way ANOVA + Tukey HSD (4 groups)
# =================================================================
endpoint <- vol %>% dplyr::group_by(mouse, group) %>%
  dplyr::filter(week == max(week)) %>% dplyr::ungroup()
cat("\n========== 8B endpoint (n per group) ==========\n"); print(table(endpoint$group))
a8 <- aov(volume ~ group, data = endpoint)
cat("\nOne-way ANOVA:\n"); print(summary(a8))
cat("\nTukey HSD:\n"); print(TukeyHSD(a8)$group)
# alternatives if skewed: kruskal.test(volume~group,endpoint); aov(log(volume)~group,endpoint)

end_summ <- endpoint %>% dplyr::group_by(group) %>%
  dplyr::summarise(mean = mean(volume), sem = sd(volume)/sqrt(dplyr::n()), n = dplyr::n(), .groups = "drop")

p8B <- ggplot(endpoint, aes(group, volume, color = group)) +
  geom_jitter(width = 0.12, size = 1.8, alpha = 0.85) +
  geom_errorbar(data = end_summ, aes(group, y = mean, ymin = mean-sem, ymax = mean+sem),
                width = 0.25, color = "black", inherit.aes = FALSE) +
  geom_point(data = end_summ, aes(group, mean), color = "black", size = 2.4, inherit.aes = FALSE) +
  scale_color_manual(values = grp_cols8) +
  labs(x = NULL, y = expression("Endpoint tumor volume (mm"^3*")")) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 25, hjust = 1))
print(p8B)
# ggsave("Fig8B_endpoint.pdf", p8B, width = 3.4, height = 3.4)

# =================================================================
# helper: pairwise LMM for a two-group contrast on a growth curve
# =================================================================
pair_lmm <- function(df, yvar, g1, g2, logscale = TRUE, label = "") {
  sub <- df %>% dplyr::filter(group %in% c(g1, g2)) %>%
    dplyr::mutate(group = droplevels(factor(group, levels = c(g1, g2))),
                  week_c = week - mean(week),
                  y = if (logscale) log(.data[[yvar]]) else .data[[yvar]])
  m <- lmer(y ~ group * week_c + (1 + week_c | mouse), data = sub, REML = TRUE,
            control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
  if (isSingular(m))
    m <- lmer(y ~ group * week_c + (1 | mouse), data = sub, REML = TRUE,
              control = lmerControl(optimizer = "bobyqa"))
  cat("\n========== LMM:", label, "(", g1, "vs", g2, ") ==========\n")
  print(anova(m, type = 3))   # 'group' = overall level diff; 'group:week_c' = divergence over time
  invisible(m)
}

# =================================================================
# S4A  TUMOR GROWTH CURVES : mean ± SEM + LMM (Am80+PD-L1 vs PD-L1)
# =================================================================
vol_summ <- vol %>% dplyr::group_by(group, week) %>%
  dplyr::summarise(mean = mean(volume), sem = sd(volume)/sqrt(dplyr::n()), n = dplyr::n(), .groups = "drop")

pair_lmm(vol, "volume", "PD-L1", "Am80 + PD-L1", logscale = TRUE, label = "S4A tumor growth")

pS4A <- ggplot(vol_summ, aes(week, mean, color = group, fill = group)) +
  geom_ribbon(aes(ymin = mean-sem, ymax = mean+sem), alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.3) +
  scale_color_manual(values = grp_cols8) + scale_fill_manual(values = grp_cols8) +
  labs(x = "Weeks", y = expression("Tumor volume (mm"^3*")"), color = NULL, fill = NULL) +
  theme_classic(base_size = 12)
print(pS4A)
# ggsave("FigS4A_growth.pdf", pS4A, width = 4.6, height = 3.2)

# =================================================================
# S4B  BODY WEIGHT CURVES : mean ± SEM + LMM (Am80+PD-L1 vs PD-L1)
# =================================================================
cat("\n=== body weight mice per group ===\n")
print(bw %>% dplyr::distinct(group, mouse) %>% dplyr::count(group))

bw_summ <- bw %>% dplyr::group_by(group, week) %>%
  dplyr::summarise(mean = mean(bw), sem = sd(bw)/sqrt(dplyr::n()), n = dplyr::n(), .groups = "drop")

pair_lmm(bw, "bw", "PD-L1", "Am80 + PD-L1", logscale = FALSE, label = "S4B body weight")

pS4B <- ggplot(bw_summ, aes(week, mean, color = group, fill = group)) +
  geom_ribbon(aes(ymin = mean-sem, ymax = mean+sem), alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.3) +
  scale_color_manual(values = grp_cols8) + scale_fill_manual(values = grp_cols8) +
  labs(x = "Weeks", y = "Body weight (g)", color = NULL, fill = NULL) +
  theme_classic(base_size = 12)
print(pS4B)
# ggsave("FigS4B_bodyweight.pdf", pS4B, width = 4.6, height = 3.2)

pS4B <- ggplot(bw_summ, aes(week, mean, color = group, fill = group)) +
  geom_ribbon(aes(ymin = mean-sem, ymax = mean+sem), alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.3) +
  scale_color_manual(values = grp_cols8) + scale_fill_manual(values = grp_cols8) +
  coord_cartesian(ylim = c(16, 24)) +                       # <- Y軸を広げる
  labs(x = "Weeks", y = "Body weight (g)", color = NULL, fill = NULL) +
  theme_classic(base_size = 12)
print(pS4B)

p_int <- "Am80+PD-L1 vs PD-L1: interaction p<0.001"   # 1.2e-5 を簡潔表記
# 厳密な値で書くなら: "interaction p=1.2×10^-5"（下記参照）

pS4A <- ggplot(vol_summ, aes(week, mean, color = group, fill = group)) +
  geom_ribbon(aes(ymin = mean-sem, ymax = mean+sem), alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.3) +
  scale_color_manual(values = grp_cols8) + scale_fill_manual(values = grp_cols8) +
  annotate("text", x = -Inf, y = Inf, label = p_int,
           hjust = -0.05, vjust = 1.5, size = 2.9) +
  labs(x = "Weeks", y = expression("Tumor volume (mm"^3*")"), color = NULL, fill = NULL) +
  theme_classic(base_size = 12)
print(pS4A)

p8B|pS4A|pS4B


com_plot<-p8B|pS4A+  theme(legend.position = "none")   |pS4B+  theme(legend.position = "none")   
ggsave("output/Fig8BS4AB.pdf", com_plot,
       width = 240,
       height = 80,
       units = "mm",
       device = cairo_pdf)
