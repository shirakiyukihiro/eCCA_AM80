# ============================================================
# Figure 6 — subcutaneous tumour, WT vs Meflin-knockout
#
# Study      : Distinct cancer-associated fibroblast lineages and an
#              Am80-based stromal-reprogramming strategy in extrahepatic
#              cholangiocarcinoma (Meflin/ISLR).
# Repository : eCCA_AM80  (https://github.com/shirakiyukihiro/eCCA_AM80)
#
# Purpose    : Tumour-growth and body-weight curves (LMM) and endpoint tumour volume (t-test) for the subcutaneous WT vs Meflin-/- experiment.
# Figure(s)  : Figure 6B/C/D.
# Input      : data/subcutaneous_implantation.xlsx.
# Output     : output/Fig6BCD.pdf.
#
# Paths have been made relative for public release: place input files
# in ./data and write outputs to ./output (create these folders, or edit
# the paths below).
# ============================================================

# ============================================================
# FIGURE 6  (subcutaneous implantation, WT vs Meflin KO)
#   6B  Tumor growth curve   -> LMM + mean±SEM plot
#   6C  Body weight curve    -> LMM + mean±SEM plot
#   6D  Endpoint tumor volume (week 8 / day 56) -> dot plot + t-test
# All error bars = SEM. All panels drawn in R (Prism retired).
# ============================================================
# Data file: subcutaneous_implantation.xlsx
#   sheet 1  'tumor (♂　8) (3)'  : tumor volumes (wide; weeks in row 1)
#   sheet 2  'BW  (♂) 計測'       : body weights (wide; weeks in row 1)
# Layout gotchas handled below:
#   - WT rows = W1..W21, KO rows = K1..K22
#   - summary rows "Ave/Max/min/SD/2SD" are interleaved -> dropped by genotype filter
#   - BW sheet has a junk column at week "9" (not body weight) -> dropped by week<=8
# ------------------------------------------------------------

library(readxl)
library(tidyverse)
library(lme4)
library(lmerTest)

file_path <- "data/subcutaneous_implantation.xlsx"

WEEK_MAX  <- 8          # last real measurement week (day 56). drops junk/empty cols.
col_wt    <- "grey40"
col_ko    <- "firebrick"

## ---- helper: read a wide sheet and return tidy long data --------
read_long <- function(path, sheet, value_name = "value", week_max = WEEK_MAX) {
  raw <- read_excel(path, sheet = sheet)
  names(raw)[1] <- "mouse"
  raw %>%
    mutate(across(-mouse, as.character)) %>%          # unify types before pivot
    pivot_longer(-mouse, names_to = "week", values_to = "v") %>%
    mutate(
      week  = suppressWarnings(as.numeric(week)),
      value = suppressWarnings(as.numeric(na_if(na_if(v, "-"), "")))
    ) %>%
    filter(!is.na(week), !is.na(value), week <= week_max) %>%
    mutate(
      genotype = case_when(
        str_starts(mouse, regex("W", ignore_case = TRUE)) ~ "WT",
        str_starts(mouse, regex("K", ignore_case = TRUE)) ~ "KO",
        TRUE ~ NA_character_                            # Ave/Max/min/SD/2SD -> NA
      )
    ) %>%
    filter(!is.na(genotype)) %>%                        # drop summary rows
    mutate(genotype = factor(genotype, levels = c("WT","KO")),
           mouse    = factor(mouse)) %>%
    transmute(mouse, genotype, week, !!value_name := value)
}

## ---- helper: fit LMM (random intercept + slope; fall back if singular) ----
fit_lmm <- function(df, yvar) {
  df <- df %>% mutate(week_c = week - mean(week),
                      y = log(.data[[yvar]]))          # log scale (var grows w/ mean)
  m_full <- lmer(y ~ genotype * week_c + (1 + week_c | mouse), data = df,
                 REML = TRUE, control = lmerControl(optimizer = "bobyqa",
                                                    optCtrl = list(maxfun = 2e5)))
  fit <- if (isSingular(m_full))
    lmer(y ~ genotype * week_c + (1 | mouse), data = df, REML = TRUE,
         control = lmerControl(optimizer = "bobyqa"))
  else m_full
  list(fit = fit, data = df)
}

## ---- helper: mean +/- SEM summary for plotting ------------------
summ_sem <- function(df, yvar) {
  df %>% group_by(genotype, week) %>%
    summarise(mean = mean(.data[[yvar]]),
              sem  = sd(.data[[yvar]]) / sqrt(n()),
              n    = n(), .groups = "drop")
}

## ---- helper: growth-curve plot with genotype P annotation -------
plot_curve <- function(summ, fit, ylab) {
  p_geno <- anova(fit, type = 3)["genotype", "Pr(>F)"]
  lab    <- paste0("Genotype effect, p = ", formatC(p_geno, format = "f", digits = 4), " (LMM)")
  ggplot(summ, aes(week, mean, color = genotype, fill = genotype)) +
    geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem), alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.8) + geom_point(size = 1.6) +
    scale_color_manual(values = c(WT = col_wt, KO = col_ko)) +
    scale_fill_manual(values  = c(WT = col_wt, KO = col_ko)) +
    annotate("text", x = min(summ$week), y = Inf, label = lab,
             hjust = 0, vjust = 1.4, size = 3.5) +
    labs(x = "Weeks", y = ylab, color = NULL, fill = NULL) +
    theme_classic(base_size = 12) +
    theme(legend.position = c(0.02, 0.9), legend.justification = c(0, 1))
}

# =================================================================
# 6B  TUMOR GROWTH
# =================================================================
tumor <- read_long(file_path, sheet = 1, value_name = "volume")
cat("\n[6B] mice per genotype:\n"); print(tapply(tumor$mouse, tumor$genotype, function(x) length(unique(x))))
cat("[6B] week range:", paste(range(tumor$week), collapse = "-"), "\n")

res_t <- fit_lmm(tumor, "volume")
cat("\n===== 6B  Tumor growth: LMM =====\n")
print(summary(res_t$fit))
print(anova(res_t$fit, type = 3))
# interaction LRT (refit both with ML, no update())
mi  <- lmer(y ~ genotype * week_c + (1 + week_c | mouse), data = res_t$data, REML = FALSE,
            control = lmerControl(optimizer = "bobyqa"))
mni <- lmer(y ~ genotype + week_c + (1 + week_c | mouse), data = res_t$data, REML = FALSE,
            control = lmerControl(optimizer = "bobyqa"))
cat("\n--- 6B interaction LRT ---\n"); print(anova(mni, mi))

p6B <- plot_curve(summ_sem(tumor, "volume"), res_t$fit,
                  expression("Tumor volume (mm"^3*")"))
print(p6B)
# ggsave("Fig6B_tumor_growth.pdf", p6B, width = 4.2, height = 3.2)

# =================================================================
# 6C  BODY WEIGHT
# =================================================================
bw <- read_long(file_path, sheet = 2, value_name = "bw")
cat("\n[6C] mice per genotype:\n"); print(tapply(bw$mouse, bw$genotype, function(x) length(unique(x))))
cat("[6C] week range:", paste(range(bw$week), collapse = "-"), "\n")  # should be 0-8, NOT 9

res_b <- fit_lmm(bw, "bw")
cat("\n===== 6C  Body weight: LMM =====\n")
print(summary(res_b$fit))
print(anova(res_b$fit, type = 3))   # expect genotype P > 0.05 (comparable)

p6C <- plot_curve(summ_sem(bw, "bw"), res_b$fit, "Body weight (g)")
print(p6C+coord_cartesian(ylim = c(22, 30)) )
# ggsave("Fig6C_body_weight.pdf", p6C, width = 4.2, height = 3.2)

# =================================================================
# 6D  ENDPOINT TUMOR VOLUME (week 8 = day 56)
# =================================================================
endpoint <- tumor %>% filter(week == WEEK_MAX)
cat("\n[6D] n per genotype at endpoint:\n"); print(table(endpoint$genotype))

# two-tailed unpaired t-test (matches current Methods/legend).
tt <- t.test(volume ~ genotype, data = endpoint, var.equal = FALSE)
cat("\n===== 6D  Endpoint t-test (Welch) =====\n"); print(tt)
# Alternatives if distribution is skewed (tumor volumes often are):
#   wilcox.test(volume ~ genotype, data = endpoint)            # Mann-Whitney
#   t.test(log(volume) ~ genotype, data = endpoint)            # t-test on log
p_end <- tt$p.value

end_summ <- endpoint %>% group_by(genotype) %>%
  summarise(mean = mean(volume), sem = sd(volume)/sqrt(n()), .groups = "drop")

p6D <- ggplot(endpoint, aes(genotype, volume, color = genotype)) +
  geom_jitter(width = 0.12, size = 1.8, alpha = 0.8) +
  geom_errorbar(data = end_summ, aes(genotype, y = mean, ymin = mean - sem, ymax = mean + sem),
                width = 0.25, inherit.aes = FALSE) +
  geom_point(data = end_summ, aes(genotype, mean), inherit.aes = FALSE, size = 2.5) +
  scale_color_manual(values = c(WT = col_wt, KO = col_ko)) +
  annotate("text", x = 1.5, y = Inf, vjust = 1.4,
           label = paste0("p = ", formatC(p_end, format = "f", digits = 4)), size = 3.5) +
  labs(x = NULL, y = expression("Tumor volume at day 56 (mm"^3*")")) +
  theme_classic(base_size = 12) + theme(legend.position = "none")
print(p6D)
# ggsave("Fig6D_endpoint.pdf", p6D, width = 2.6, height = 3.2)

# =================================================================
# DIAGNOSTICS (check residuals for the two LMMs)
# =================================================================
op <- par(mfrow = c(2,2))
plot(fitted(res_t$fit), resid(res_t$fit), main="6B resid"); abline(h=0, lty=2)
qqnorm(resid(res_t$fit), main="6B QQ"); qqline(resid(res_t$fit))
plot(fitted(res_b$fit), resid(res_b$fit), main="6C resid"); abline(h=0, lty=2)
qqnorm(resid(res_b$fit), main="6C QQ"); qqline(resid(res_b$fit))
par(op)

com_plot<-p6B|(p6C+coord_cartesian(ylim = c(22, 30)) )|p6D
ggsave("output/Fig6BCD.pdf", com_plot,
       width = 240,
       height = 70,
       units = "mm",
       device = cairo_pdf)
