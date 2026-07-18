# ============================================================
# Subcutaneous tumour IHC quantification — WT vs Meflin-knockout
#
# Study      : Distinct cancer-associated fibroblast lineages and an
#              Am80-based stromal-reprogramming strategy in extrahepatic
#              cholangiocarcinoma (Meflin/ISLR).
# Repository : eCCA_AM80  (https://github.com/shirakiyukihiro/eCCA_AM80)
#
# Purpose    : Read subcutaneous-tumour immunostaining quantifications and compare wild-type versus Meflin-knockout (LMM-based analysis and plots). Shares the input workbook with the Am80 IHC scripts.
# Figure(s)  : Figure 6E-H (subcutaneous IHC box plots, WT vs Meflin-/-).
# Input      : data/subcutaneous_IHC_data.xlsx.
# Output     : output/2026.06.13_KO_quantification.pdf.
#
# Paths have been made relative for public release: place input files
# in ./data and write outputs to ./output (create these folders, or edit
# the paths below).
# ============================================================

list.files("data")
file_path <- "data/subcutaneous_IHC_data.xlsx"

#file_path <- "data/subcutaneous_IHC_WT_KO.xlsx"
sheets<-excel_sheets(file_path)
# Sheet1を読み込み
df <- read_excel(file_path, sheet = sheets[1],skip = 1)
colnames <- read_excel(file_path, sheet = sheets[1],n_max = 1)

# 確認
head(df)
colnames(colnames)

library(readxl)
library(tidyverse)
library(lme4)
library(lmerTest)
library(ggbeeswarm)
library(patchwork)


##################################################
# Convert to long format
##################################################

marker_info <- tribble(
  ~Marker,   ~WT_col,   ~KO_col,
  "Ki67",    "WT...3",  "KO...4",
  "CD31",    "WT...6",  "KO...7",
  "SMA",     "WT...9",  "KO...10",
  "MT",      "WT...12", "KO...13",
  "Pdgfrb",  "WT...15", "KO...16",
  "CD11b",   "WT...18", "KO...19"
)


count_table <- map_dfr(
  1:nrow(marker_info),
  function(i){
    
    tibble(
      Marker = marker_info$Marker[i],
      
      WT_n = sum(
        !is.na(df[[marker_info$WT_col[i]]])
      ),
      
      KO_n = sum(
        !is.na(df[[marker_info$KO_col[i]]])
      )
    )
  }
)

count_table

##################################################
# Mouse structure
##################################################

mouse_structure <- list(
  
  Ki67 = list(
    WT = rep(5, 21),
    KO = rep(5, 22)
  ),
  
  CD31 = list(
    WT = rep(5, 21),
    KO = rep(5, 22)
  ),
  
  SMA = list(
    WT = c(
      rep(5, 3),
      4,
      rep(5, 7),
      4,
      rep(5, 9)
    ),
    KO = rep(5, 22)
  ),
  
  MT = list(
    WT = c(
      rep(5, 7),
      3,
      rep(5, 3),
      2,
      rep(5, 9)
    ),
    KO = rep(5, 22)
  ),
  
  Pdgfrb = list(
    WT = c(
      rep(5, 3),
      3,
      rep(5, 3),
      3,
      rep(5, 13)
    ),
    KO = rep(5, 22)
  ),
  
  CD11b = list(
    WT = c(
      rep(5, 3),
      4,
      rep(5, 17)
    ),
    KO = c(
      rep(5, 3),
      4,
      rep(5, 18)
    )
  )
)

##################################################
# Helper function
##################################################

make_marker_df <- function(marker, wt_col, ko_col){
  
  wt_values <- df[[wt_col]]
  wt_values <- wt_values[!is.na(wt_values)]
  
  ko_values <- df[[ko_col]]
  ko_values <- ko_values[!is.na(ko_values)]
  
  wt_counts <- mouse_structure[[marker]]$WT
  ko_counts <- mouse_structure[[marker]]$KO
  
  wt_df <- tibble(
    Marker = marker,
    Group = "WT",
    MouseID = rep(
      paste0("WT_", seq_along(wt_counts)),
      times = wt_counts
    ),
    Value = wt_values
  )
  
  ko_df <- tibble(
    Marker = marker,
    Group = "KO",
    MouseID = rep(
      paste0("KO_", seq_along(ko_counts)),
      times = ko_counts
    ),
    Value = ko_values
  )
  
  bind_rows(wt_df, ko_df)
}

##################################################
# Long format with correct MouseID
##################################################

df_long <- bind_rows(
  
  make_marker_df(
    "Ki67",
    "WT...3",
    "KO...4"
  ),
  
  make_marker_df(
    "CD31",
    "WT...6",
    "KO...7"
  ),
  
  make_marker_df(
    "SMA",
    "WT...9",
    "KO...10"
  ),
  
  make_marker_df(
    "MT",
    "WT...12",
    "KO...13"
  ),
  
  make_marker_df(
    "Pdgfrb",
    "WT...15",
    "KO...16"
  ),
  
  make_marker_df(
    "CD11b",
    "WT...18",
    "KO...19"
  )
)

##################################################
# Sanity check
##################################################

df_long %>%
  dplyr::count(Marker, Group, MouseID) %>%
  arrange(Marker, Group, MouseID)
##################################################
# Mouse mean for plotting
##################################################

mouse_mean <- df_long %>%
  group_by(
    Marker,
    Group,
    MouseID
  ) %>%
  summarise(
    Value = mean(Value),
    .groups = "drop"
  )

##################################################
# Mixed model statistics
##################################################

library(lme4)
library(lmerTest)

stats <- purrr::map_dfr(
  unique(df_long$Marker),
  function(marker){
    
    dat <- df_long %>%
      filter(Marker == marker)
    
    fit <- lmer(
      Value ~ Group + (1|MouseID),
      data = dat
    )
    
    pval <- anova(fit)$`Pr(>F)`[1]
    
    tibble(
      Marker = marker,
      p = pval
    )
  }
)

print(stats)

##################################################
# Colors
##################################################

cols <- c(
  WT = "#18499e",
  KO = "#e7211a"
)

##################################################
# Theme
##################################################

theme_my <- theme_classic(base_size = 10) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text = element_text(
      color = "black",
      size = 9
    ),
    axis.line = element_line(
      color = "black"
    ),
    plot.title = element_text(
      hjust = 0.5,
      size = 10
    )
  )

##################################################
# Plot function
##################################################

library(ggbeeswarm)

plot_marker <- function(marker_name){
  
  dat <- mouse_mean %>%
    filter(Marker == marker_name)
  
  pval <- stats %>%
    filter(Marker == marker_name) %>%
    pull(p)
  
  ymax <- max(dat$Value, na.rm = TRUE)
  
  ggplot(
    dat,
    aes(
      x = Group,
      y = Value
    )
  ) +
    
    geom_boxplot(
      aes(
        fill = Group
      ),
      width = 0.30,
      alpha = 0.35,
      color = "black",
      linewidth = 0.8,
      outlier.shape = NA,
      staplewidth = 0.5
    ) +
    
    geom_beeswarm(
      color = "black",
      size = 1.6,
      alpha = 0.4,
      priority = "density",
      cex = 5
    ) +
    
    scale_fill_manual(
      values = cols
    ) +
    
    annotate(
      "text",
      x = 1.5,
      y = ymax * 1.08,
      label = paste0(
        "p = ",
        signif(pval, 3)
      ),
      size = 3
    ) +
    
    coord_cartesian(
      ylim = c(
        min(dat$Value, na.rm = TRUE),
        ymax * 1.15
      )
    ) +
    
    labs(
      title = marker_name
    ) +
    
    theme_my
}

##################################################
# Create plots
##################################################

p1 <- plot_marker("Ki67")
p2 <- plot_marker("CD31")
p3 <- plot_marker("SMA")
p4 <- plot_marker("MT")
p5 <- plot_marker("Pdgfrb")
p6 <- plot_marker("CD11b")

##################################################
# Combine
##################################################

library(patchwork)

combined_plot <-
  (p4 +p5  + p3) /
  (p1 + p2 + p6)

combined_plot

##################################################
# Save
##################################################

ggsave(
  filename = "output/2026.06.13_KO_quantification.pdf",
  plot = combined_plot,
  width = 180,
  height = 120,
  units = "mm",
  device = cairo_pdf
)

