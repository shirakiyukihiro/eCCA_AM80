# ============================================================
# 3D tissue clearing — lineage-traced cell morphometry (elongation / sphericity)
#
# Study      : Distinct cancer-associated fibroblast lineages and an
#              Am80-based stromal-reprogramming strategy in extrahepatic
#              cholangiocarcinoma (Meflin/ISLR).
# Repository : eCCA_AM80  (https://github.com/shirakiyukihiro/eCCA_AM80)
#
# Purpose    : Read Imaris-derived per-cell morphometry (elongation index and sphericity) of lineage-traced cells and compare stromal vs biliary compartments.
# Figure(s)  : 3D morphometry panel (elongation index / sphericity; verify figure number against final legends).
# Input      : data/clearing_morphometry.xlsx.
# Output     : output/2026.06.03_Elongation_Sphericity.pdf.
#
# Paths have been made relative for public release: place input files
# in ./data and write outputs to ./output (create these folders, or edit
# the paths below).
# ============================================================

list.files("data")
library(readxl)

file_path <- "data/clearing_morphometry.xlsx"

# Sheet1を読み込み
df <- read_excel(file_path, sheet = 1)

# 確認
head(df)
colnames(df)
library(tidyverse)
library(ggbeeswarm)
library(tidyverse)
library(ggbeeswarm)
library(patchwork)

# グループ名変更
df <- df %>%
  mutate(
    area = factor(
      area,
      levels = c("B", "S"),
      labels = c("Biliary", "Stromal")
    )
  )

cols <- c(
  "Biliary" = "#18499e",
  "Stromal" = "#e7211a"
)

theme_my <- theme_classic(base_size = 3) +
  theme(
    legend.position = "none",
    axis.title = element_blank(),
    axis.line = element_line(
      color = "black",
      linewidth = 1
    ),
    axis.text = element_text(
      size = 3,
      color = "black"
    )
  )

##################################################
# Elongation Index
##################################################

p1 <- ggplot()+
  
  geom_beeswarm(
    data = subset(df, area == "Biliary"),
    aes(x = area, y = `Elongation Index`),
    shape = 16,
    size = 0.8,
    color = cols["Biliary"],
    alpha = 0.1,
    priority = "density",
    cex = 1
  ) +geom_beeswarm(
    data = subset(df, area == "Stromal"),
    aes(x = area, y = `Elongation Index`),
    shape = 15,
    size = 1,
    color = cols["Stromal"],
    alpha = 0.3,
    priority = "density",
    cex = 1
  )   +
  
  geom_boxplot(
    data = df,
    aes(x = area, y = `Elongation Index`),
    width = 0.3,
    fill = NA,
    color = "black",
    outlier.shape = NA,
    linewidth = 1,
    staplewidth = 0.5
  ) +
  
  ylab("Elongation Index") +
  theme_my

p1+
  scale_y_log10()
##################################################
# Sphericity
##################################################

p2 <- ggplot()+
  
  geom_beeswarm(
    data = subset(df, area == "Biliary"),
    aes(x = area, y = Sphericity),
    shape = 16,
    size = 0.8,
    color = cols["Biliary"],
    alpha = 0.1,
    priority = "density",
    cex = 1
  ) +
  geom_beeswarm(
    data = subset(df, area == "Stromal"),
    aes(x = area, y = Sphericity),
    shape = 15,
    size = 1,
    color = cols["Stromal"],
    alpha = 0.3,
    priority = "density",
    cex = 1
  ) +
  
  geom_boxplot(
    data = df,
    aes(x = area, y = Sphericity),
    width = 0.3,
    fill = NA,
    color = "black",
    outlier.shape = NA,
    linewidth = 1,
    staplewidth = 0.5
  ) +
  
  ylab("Sphericity") +
  theme_my

##################################################
# combine
##################################################

p1 
p2

library(patchwork)

p1 <- p1 +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    )
  )+  scale_y_log10()

p2 <- p2 +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    )
  )

library(patchwork)
combined_plot <- p1 | p2 


ggsave(
  filename = "output/2026.06.03_Elongation_Sphericity.pdf",
  plot = combined_plot,
  width = 90,
  height = 50,
  units = "mm",
  device = cairo_pdf
)

