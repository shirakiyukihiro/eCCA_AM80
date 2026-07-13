# ============================================================
# Figure 7 — Am80 + chemotherapy (gemcitabine / cisplatin)
#
# Study      : Distinct cancer-associated fibroblast lineages and an
#              Am80-based stromal-reprogramming strategy in extrahepatic
#              cholangiocarcinoma (Meflin/ISLR).
# Repository : eCCA_AM80  (https://github.com/shirakiyukihiro/eCCA_AM80)
#
# Purpose    : Endpoint tumour volume (one-way ANOVA + Tukey HSD within genotype) and tumour-growth / body-weight curves for the Am80 + chemotherapy experiment.
# Figure(s)  : Figure 7B/C/D.
# Input      : data/Am80_chemotherapy.xlsx.
# Output     : output/Fig7BCD.pdf.
#
# Paths have been made relative for public release: place input files
# in ./data and write outputs to ./output (create these folders, or edit
# the paths below).
# ============================================================

# ============================================================
# FIGURE 7  (Am80 + chemotherapy; gemcitabine/cisplatin)
#   7B  Endpoint tumor volume  -> one-way ANOVA + Tukey HSD, within each genotype
#   7C  Tumor growth curves    -> mean ± SEM (descriptive)
#   7D  Body weight curves      -> mean ± SEM (descriptive)
# All error bars = SEM. WT = 4 groups, KO = 3 groups.
# ============================================================
# Data file: Am80_chemotherapy.xlsx
#   sheet 'WT'     : WT tumor volumes (4 groups)
#   sheet 'KO (2)' : KO tumor volumes (3 groups)
#   sheet 'WT-BW'  : WT body weights
#   sheet 'KO-BW'  : KO body weights
#
# Treatment label (col 1)  ->  Figure group:
#   DMSO/PBS -> Control        Am80/PBS -> Am80
#   DMSO/CT  -> Chemo          Am80/CT  -> Am80 + Chemo      (CT = gem/cis)
#
# Sheet quirks handled:
#   - week headers carry annotations e.g. "4(Am80)","5(CT1)" -> stripped to numeric
#   - WT weeks 4-9, KO weeks 0-5 -> re-based to "weeks since first measurement" (0-5)
#   - NO mouse-ID column -> IDs assigned by row within sheet
#   - summary rows (group mean + SD) sit BELOW the data with the SAME group labels
#     -> data block detected as the first contiguous appearance of each group
# ------------------------------------------------------------

library(tidyverse)
library(conflicted)
for (f in c("select","filter","mutate","rename","count","distinct",
            "summarise","summarize","group_by","arrange","transmute","slice"))
  conflict_prefer(f, "dplyr")
conflict_prefer("pivot_longer", "tidyr")

file_path <- "data/Am80_chemotherapy.xlsx"    

valid_groups <- c("Am80/CT","Am80/PBS","DMSO/CT","DMSO/PBS")
grp_map <- c("DMSO/PBS"="Control","Am80/PBS"="Am80","DMSO/CT"="Chemo","Am80/CT"="Am80 + Chemo")
grp_levels <- c("Control","Am80","Chemo","Am80 + Chemo")
grp_cols   <- c("Control"="grey55","Am80"="#377EB8","Chemo"="#FF7F00","Am80 + Chemo"="#E41A1C")

## ---- read one sheet, extract the data block, return tidy long ----
read_block <- function(path, sheet, value_name = "value") {
  raw <- read_excel(path, sheet = sheet, col_names = FALSE)
  # row 1 = week header (col 1 blank). Parse week numbers from row 1.
  hdr <- as.character(unlist(raw[1, ]))
  wk  <- suppressWarnings(as.numeric(gsub("\\(.*\\)", "", hdr)))  # "4(Am80)"->4
  body <- raw[-1, ]
  grp_raw <- as.character(body[[1]])
  
  # --- detect the data block: first contiguous appearance of each group ---
  keep <- logical(length(grp_raw)); seen <- character(0); cur <- NA_character_
  for (k in seq_along(grp_raw)) {
    g <- grp_raw[k]
    if (is.na(g) || !(g %in% valid_groups)) break
    if (!is.na(cur) && g != cur && g %in% seen) break
    if (is.na(cur) || g != cur) { if (!is.na(cur)) seen <- c(seen, cur); cur <- g }
    keep[k] <- TRUE
  }
  dat <- body[keep, , drop = FALSE]
  grp <- grp_raw[keep]
  
  # column 1 = group label; columns 2..n = week measurements
  meas <- dat[, -1, drop = FALSE]                 # measurement columns only
  meas <- dplyr::mutate(meas, dplyr::across(dplyr::everything(), as.character))
  wk_cols <- wk[-1]                               # corresponding week numbers
  colnames(meas) <- as.character(wk_cols)
  meas$.mouse <- paste0(sheet, "_", seq_len(nrow(meas)))
  meas$.gid   <- grp
  
  long <- meas %>%
    tidyr::pivot_longer(cols = -c(.mouse, .gid),
                        names_to = "week", values_to = "v") %>%
    dplyr::mutate(week  = suppressWarnings(as.numeric(week)),
                  value = suppressWarnings(as.numeric(v))) %>%
    dplyr::filter(!is.na(week), !is.na(value)) %>%
    dplyr::mutate(week_rel = week - min(week, na.rm = TRUE),
                  group = factor(unname(grp_map[.gid]), levels = grp_levels),
                  mouse = factor(.mouse)) %>%
    dplyr::transmute(mouse, group, week, week_rel, !!value_name := value)
  long
}

# =================================================================
# Load tumor volume (WT, KO) and body weight (WT-BW, KO-BW)
# =================================================================
wt   <- read_block(file_path, "WT",     "volume") %>% mutate(genotype = "WT")
ko   <- read_block(file_path, "KO (2)", "volume") %>% mutate(genotype = "KO")
wtbw <- read_block(file_path, "WT-BW",  "bw")     %>% mutate(genotype = "WT")
kobw <- read_block(file_path, "KO-BW",  "bw")     %>% mutate(genotype = "KO")

vol <- bind_rows(wt, ko) %>% mutate(genotype = factor(genotype, levels = c("WT","KO")))
bw  <- bind_rows(wtbw, kobw) %>% mutate(genotype = factor(genotype, levels = c("WT","KO")))

# ---- SANITY CHECKS: verify group sizes vs the legend ----
cat("\n=== mice per group (volume) — expect WT 5/5/5/5, KO Control5 Chemo4 Am80+Chemo6 ===\n")
print(vol %>% distinct(genotype, group, mouse) %>% count(genotype, group))
cat("\n=== week_rel range per genotype (expect 0-5) ===\n")
print(vol %>% group_by(genotype) %>% summarise(min = min(week_rel), max = max(week_rel)))

# =================================================================
# 7B  ENDPOINT TUMOR VOLUME : one-way ANOVA + Tukey, within genotype
# =================================================================
endpoint <- vol %>% group_by(genotype, mouse, group) %>%
  filter(week_rel == max(week_rel)) %>% ungroup()      # last timepoint per mouse

run_anova <- function(df, label) {
  df <- droplevels(df)
  cat("\n========== 7B endpoint:", label, "==========\n")
  cat("n per group:\n"); print(table(df$group))
  a <- aov(volume ~ group, data = df)
  cat("\nOne-way ANOVA:\n"); print(summary(a))
  cat("\nTukey HSD:\n"); print(TukeyHSD(a)$group)
  # alternatives if assumptions are violated:
  #   kruskal.test(volume ~ group, df); aov(log(volume) ~ group, df)
  invisible(a)
}
run_anova(filter(endpoint, genotype == "WT"), "WT (4 groups)")
run_anova(filter(endpoint, genotype == "KO"), "KO (3 groups)")

end_summ <- endpoint %>% group_by(genotype, group) %>%
  summarise(mean = mean(volume), sem = sd(volume)/sqrt(n()), n = n(), .groups = "drop")

p7B <- ggplot(endpoint, aes(group, volume, color = group)) +
  geom_jitter(width = 0.12, size = 1.6, alpha = 0.8) +
  geom_errorbar(data = end_summ, aes(group, y = mean, ymin = mean-sem, ymax = mean+sem),
                width = 0.25, color = "black", inherit.aes = FALSE) +
  geom_point(data = end_summ, aes(group, mean), color = "black", size = 2.3, inherit.aes = FALSE) +
  facet_wrap(~ genotype, scales = "free_x") +
  scale_color_manual(values = grp_cols) +
  labs(x = NULL, y = expression("Endpoint tumor volume (mm"^3*")")) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
print(p7B)
# Add significance brackets from Tukey with ggsignif if desired:
#   library(ggsignif); + geom_signif(comparisons=list(c("Chemo","Am80 + Chemo")), ...)
# ggsave("Fig7B_endpoint.pdf", p7B, width = 6, height = 3.4)

# =================================================================
# 7C  TUMOR GROWTH CURVES : mean ± SEM (descriptive)
# =================================================================
vol_summ <- vol %>% group_by(genotype, group, week_rel) %>%
  summarise(mean = mean(volume), sem = sd(volume)/sqrt(n()), n = n(), .groups = "drop")

p7C <- ggplot(vol_summ, aes(week_rel, mean, color = group, fill = group)) +
  geom_ribbon(aes(ymin = mean-sem, ymax = mean+sem), alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.3) +
  facet_wrap(~ genotype) +
  scale_color_manual(values = grp_cols) + scale_fill_manual(values = grp_cols) +
  labs(x = "Weeks after treatment start", y = expression("Tumor volume (mm"^3*")"),
       color = NULL, fill = NULL) +
  theme_classic(base_size = 12)
print(p7C)
# ggsave("Fig7C_growth.pdf", p7C, width = 7, height = 3.4)

# =================================================================
# 7D  BODY WEIGHT CURVES : mean ± SEM (descriptive)
# =================================================================
cat("\n=== mice per group (body weight) ===\n")
print(bw %>% distinct(genotype, group, mouse) %>% count(genotype, group))

bw_summ <- bw %>% group_by(genotype, group, week_rel) %>%
  summarise(mean = mean(bw), sem = sd(bw)/sqrt(n()), n = n(), .groups = "drop")

p7D <- ggplot(bw_summ, aes(week_rel, mean, color = group, fill = group)) +
  geom_ribbon(aes(ymin = mean-sem, ymax = mean+sem), alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.3) +
  facet_wrap(~ genotype) +
  scale_color_manual(values = grp_cols) + scale_fill_manual(values = grp_cols) +
  labs(x = "Weeks after treatment start", y = "Body weight (g)", color = NULL, fill = NULL) +
  theme_classic(base_size = 12)
print(p7D)
# ggsave("Fig7D_bodyweight.pdf", p7D, width = 7, height = 3.4)

conflicts_prefer(lmerTest::lmer)
# --- pair_lmm 関数（figure8 スクリプトと同じもの。未定義ならコピー）---
pair_lmm <- function(df, yvar, g1, g2, logscale = TRUE, label = "") {
  sub <- df %>% dplyr::filter(group %in% c(g1, g2)) %>%
    dplyr::mutate(group = droplevels(factor(group, levels = c(g1, g2))),
                  week_c = week_rel - mean(week_rel),
                  y = if (logscale) log(.data[[yvar]]) else .data[[yvar]])
  m <- lmer(y ~ group * week_c + (1 + week_c | mouse), data = sub, REML = TRUE,
            control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
  if (isSingular(m)) m <- lmer(y ~ group * week_c + (1 | mouse), data = sub, REML = TRUE,
                               control = lmerControl(optimizer = "bobyqa"))
  cat("\n===== LMM:", label, "(", g1, "vs", g2, ") =====\n"); print(anova(m, type = 3)); invisible(m)
}

# 7C 腫瘍増殖：各遺伝子型内で Am80+Chemo vs Chemo（主要対比）
pair_lmm(filter(vol, genotype=="WT"), "volume", "Chemo", "Am80 + Chemo", TRUE, "7C tumor WT")
pair_lmm(filter(vol, genotype=="KO"), "volume", "Chemo", "Am80 + Chemo", TRUE, "7C tumor KO")

# 7D 体重：同対比（安全性、有意差なしを確認）
pair_lmm(filter(bw, genotype=="WT"), "bw", "Chemo", "Am80 + Chemo", FALSE, "7D BW WT")
pair_lmm(filter(bw, genotype=="KO"), "bw", "Chemo", "Am80 + Chemo", FALSE, "7D BW KO")


p7D <- ggplot(bw_summ, aes(week_rel, mean, color = group, fill = group)) +
  geom_ribbon(aes(ymin = mean-sem, ymax = mean+sem), alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.3) +
  facet_wrap(~ genotype) +
  scale_color_manual(values = grp_cols) + scale_fill_manual(values = grp_cols) +
  coord_cartesian(ylim = c(16, 24)) +                       # <- Y軸を広げる
  labs(x = "Weeks after treatment start", y = "Body weight (g)", color = NULL, fill = NULL) +
  theme_classic(base_size = 12)
print(p7D)

library(ggsignif)
# 主要比較と Tukey p 値（先の結果から）を手動指定
sig7B <- data.frame(
  genotype = factor(c("WT","WT","KO"), levels = c("WT","KO")),
  group    = "Chemo",  # ダミー（facet 認識用）
  start    = c("Control","Chemo","Chemo"),
  end      = c("Am80 + Chemo","Am80 + Chemo","Am80 + Chemo"),
  y        = c(780, 700, 620),          # ブラケットの高さ（データに応じ調整）
  label    = c("p<0.001","p=0.046","p=0.908")
)
p7B <- p7B +
  ggsignif::geom_signif(
    data = sig7B,
    aes(xmin = start, xmax = end, annotations = label, y_position = y),
    manual = TRUE, textsize = 3, tip_length = 0.01, inherit.aes = FALSE)
print(p7B)


# 主要対比（Am80+Chemo vs Chemo）の交互作用 p 値をファセット別に注記
note7C <- data.frame(
  genotype = factor(c("WT","KO"), levels = c("WT","KO")),
  label    = c("Am80+Chemo vs Chemo\ninteraction p=0.0065",
               "Am80+Chemo vs Chemo\ninteraction p=0.824")
)

p7C <- ggplot(vol_summ, aes(week_rel, mean, color = group, fill = group)) +
  geom_ribbon(aes(ymin = mean-sem, ymax = mean+sem), alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.3) +
  facet_wrap(~ genotype) +
  scale_color_manual(values = grp_cols) + scale_fill_manual(values = grp_cols) +
  geom_text(data = note7C, aes(x = -Inf, y = Inf, label = label),
            hjust = -0.05, vjust = 1.2, size = 2.9, inherit.aes = FALSE, lineheight = 0.9) +
  labs(x = "Weeks after treatment start", y = expression("Tumor volume (mm"^3*")"),
       color = NULL, fill = NULL) +
  theme_classic(base_size = 12)
print(p7C)


com_plot<-p7B|p7C +  theme(legend.position = "none")   |p7D +  theme(legend.position = "none")   
ggsave("output/Fig7BCD.pdf", com_plot,
       width = 300,
       height = 80,
       units = "mm",
       device = cairo_pdf)
