# ============================================================
# In situ hybridisation quantification — Meflin (Islr) and Thy1 (normal vs BDL)
#
# Study      : Distinct cancer-associated fibroblast lineages and an
#              Am80-based stromal-reprogramming strategy in extrahepatic
#              cholangiocarcinoma (Meflin/ISLR).
# Repository : eCCA_AM80  (https://github.com/shirakiyukihiro/eCCA_AM80)
#
# Purpose    : Read ISH signal-intensity quantifications for Islr (Meflin) and Thy1 in the normal versus bile-duct-ligation (cholangitis) model and generate the comparison plots.
# Figure(s)  : ISH Islr/Thy1 panels (Fig. 2B/D; verify against final legends).
# Input      : data/ISH_Islr_Thy1.xlsx.
# Output     : output/2026.06.22_Fig2bd_Islr_Thy1.pdf.
#
# Paths have been made relative for public release: place input files
# in ./data and write outputs to ./output (create these folders, or edit
# the paths below).
# ============================================================

library(readxl)
library(tidyverse)
library(ggbeeswarm)
library(patchwork)

##################################################
# Read data
##################################################

file_path <- "data/ISH_Islr_Thy1.xlsx"

##################################################
# Biliary
##################################################

df <- read_excel(file_path, sheet = "biliary")

df_biliary <- bind_rows(
  df %>%
    dplyr::select(
      Group = ISLR,
      Value = `...2`
    ) %>%
    mutate(Marker = "Islr"),
  
  df %>%
    dplyr::select(
      Group = Thy1,
      Value = `...6`
    ) %>%
    mutate(Marker = "Thy1")
) %>%
  relocate(Marker) %>%
  filter(!is.na(Value)) %>%
  mutate(Region = "Biliary")

##################################################
# Stromal
##################################################

df <- read_excel(file_path, sheet = "stromal")

df_stromal <- bind_rows(
  df %>%
    dplyr::select(
      Group = Islr,
      Value = `...2`
    ) %>%
    mutate(Marker = "Islr"),
  
  df %>%
    dplyr::select(
      Group = Thy1,
      Value = `...7`
    ) %>%
    mutate(Marker = "Thy1")
) %>%
  relocate(Marker) %>%
  filter(!is.na(Value)) %>%
  mutate(Region = "Stromal")

##################################################
# Merge
##################################################

df_long <- bind_rows(
  df_biliary,
  df_stromal
)%>%mutate(Group=factor(Group,levels = c("normal","CBDL")))


## ============================ MOUSE ID ===============================
## Per-mouse field counts (IN ROW ORDER) for each Region|Marker|Group.
mouse_counts <- list(
  "Biliary|Islr|CBDL"   = c(6, 6, 5),
  "Biliary|Islr|normal" = c(5, 4, 4, 4, 4),
  "Biliary|Thy1|CBDL"   = c(5, 5, 5),
  "Biliary|Thy1|normal" = c(2, 2, 2, 2, 2),
  "Stromal|Islr|CBDL"   = c(5, 5, 5),
  "Stromal|Islr|normal" = c(5, 4, 4, 4, 4),
  "Stromal|Thy1|CBDL"   = c(5, 5, 5),
  "Stromal|Thy1|normal" = c(2, 2, 2, 1, 1)
)

assign_mouse <- function(region, marker, group, n) {
  key <- paste(region, marker, group, sep = "|")
  cnt <- mouse_counts[[key]]
  if (is.null(cnt)) stop("no counts for ", key)
  if (sum(cnt) != n)
    warning(sprintf("%s: counts sum=%d but data n=%d -> last mouse adjusted; verify.",
                    key, sum(cnt), n))
  ids <- rep(seq_along(cnt), cnt)
  if (length(ids) > n) ids <- ids[seq_len(n)]                          # truncate
  if (length(ids) < n) ids <- c(ids, rep(length(cnt), n - length(ids))) # pad
  paste0(group, "_m", ids)
}


df_long <- df_long %>%
  group_by(Region, Marker, Group) %>%
  mutate(Mouse = assign_mouse(Region[1], Marker[1],
                              as.character(Group[1]), dplyr::n())) %>%
  ungroup()


## sanity: mice per cell
base::print(as.data.frame(
  df_long %>% distinct(Region, Marker, Group, Mouse) %>%
    dplyr::count(Region, Marker, Group, name = "n_mice")
))

## ----------------------------- model + per-mouse means --------------
ylab_for <- function(mk) bquote(italic(.(mk))^"+" ~ "signal (per field)")  # <-edit label if needed

panel <- function(rg, mk) {
  d <- df_long %>% filter(Region == rg, Marker == mk)
  
  ## mixed model on FIELD-level data; p for Group (CBDL vs normal)
  m <- lmerTest::lmer(Value ~ Group + (1 | Mouse), data = d)
  p <- summary(m)$coefficients["GroupCBDL", "Pr(>|t|)"]
  cat(sprintf("%-8s %-5s  n(mouse): normal=%d CBDL=%d  p(lmer)=%.4g  %s\n",
              rg, mk, n_distinct(d$Mouse[d$Group == "normal"]),
              n_distinct(d$Mouse[d$Group == "CBDL"]), p,
              if (isSingular(m)) "[singular fit]" else ""))
  
  ## dots = per-mouse means
  pm <- d %>% group_by(Group, Mouse) %>%
    summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop")
  
  plot_box(pm, "Group", "Value", "Group",
           ylab    = ylab_for(mk),
           pval    = p,
           cols    = c(normal = "#4d4d4d", CBDL = "#850000"),
           xlabels = c(normal = "Normal", CBDL = "CBDL"),
           p_x     = 1.5)
}

## order: Biliary-Islr, Stromal-Islr, Biliary-Thy1, Stromal-Thy1
combined <-
  panel("Biliary", "Islr") + panel("Stromal", "Islr") +
  panel("Biliary", "Thy1") + panel("Stromal", "Thy1") +
  plot_layout(nrow = 1)
combined


ggsave(
  filename = "output/2026.06.22_Fig2bd_Islr_Thy1.pdf",
  plot = combined,
  width = 180,
  height = 60,
  units = "mm",
  device = cairo_pdf
)
