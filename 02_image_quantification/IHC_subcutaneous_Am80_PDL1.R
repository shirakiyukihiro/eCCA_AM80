# ============================================================
# Subcutaneous tumour IHC quantification — Am80 + anti-PD-L1
#
# Study      : Distinct cancer-associated fibroblast lineages and an
#              Am80-based stromal-reprogramming strategy in extrahepatic
#              cholangiocarcinoma (Meflin/ISLR).
# Repository : eCCA_AM80  (https://github.com/shirakiyukihiro/eCCA_AM80)
#
# Purpose    : Read subcutaneous-tumour immunostaining quantifications for the Am80 + anti-PD-L1 arms (LMM-based analysis and plots). Shares the input workbook with the WT-vs-KO and Am80/chemo IHC scripts.
# Figure(s)  : Figure 8C-L (subcutaneous IHC box plots, Am80 + anti-PD-L1).
# Input      : data/subcutaneous_IHC_data.xlsx.
# Output     : output/2026.06.13_AM80_PDL1_analysis.pdf.
#
# Paths have been made relative for public release: place input files
# in ./data and write outputs to ./output (create these folders, or edit
# the paths below).
# ============================================================

library(readxl)
library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggbeeswarm)
library(patchwork)

file_path <- "data/subcutaneous_IHC_data.xlsx"

sheets <- excel_sheets(file_path)

df <- read_excel(
  file_path,
  sheet = sheets[3],
  skip = 1
)
colnames <- read_excel(file_path, sheet = sheets[3],n_max = 1)

##################################################
# Marker information
##################################################

marker_info <- tibble(
  
  Marker = c(
    "Ki67",
    "CD31",
    "SMA",
    "MT",
    "CD3",
    "CD4",
    "CD8",
    "CD11b",
    "CD86",
    "CD163"
  ),
  
  AP = c(
    "AP...1",
    "AP...6",
    "AP...11",
    "AP...16",
    "AP...21",
    "AP...26",
    "AP...31",
    "AP...36",
    "AP...41",
    "AP...46"
  ),
  
  DP = c(
    "DP...2",
    "DP...7",
    "DP...12",
    "DP...17",
    "DP...22",
    "DP...27",
    "DP...32",
    "DP...37",
    "DP...42",
    "DP...47"
  ),
  
  AI = c(
    "AI...3",
    "AI...8",
    "AI...13",
    "AI...18",
    "AI...23",
    "AI...28",
    "AI...33",
    "AI...38",
    "AI...43",
    "AI...48"
  ),
  
  DI = c(
    "DI...4",
    "DI...9",
    "DI...14",
    "DI...19",
    "DI...24",
    "DI...29",
    "DI...34",
    "DI...39",
    "DI...44",
    "DI...49"
  )
)

##################################################
# Helper
##################################################

make_group_df <- function(
    values,
    marker,
    group){
  
  values <- values[!is.na(values)]
  
  n_mouse <- ceiling(length(values)/5)
  
  mouse_id <- rep(
    paste0(group,"_",seq_len(n_mouse)),
    each = 5
  )
  
  tibble(
    Marker = marker,
    Group = group,
    MouseID = mouse_id[seq_along(values)],
    Value = values
  )
}

##################################################
# Long format
##################################################

df_long <- map_dfr(
  seq_len(nrow(marker_info)),
  function(i){
    
    bind_rows(
      
      make_group_df(
        df[[marker_info$AP[i]]],
        marker_info$Marker[i],
        "AP"
      ),
      
      make_group_df(
        df[[marker_info$DP[i]]],
        marker_info$Marker[i],
        "DP"
      ),
      
      make_group_df(
        df[[marker_info$AI[i]]],
        marker_info$Marker[i],
        "AI"
      ),
      
      make_group_df(
        df[[marker_info$DI[i]]],
        marker_info$Marker[i],
        "DI"
      )
      
    )
    
  }
)

##################################################
# Mouse mean for plotting
##################################################

##################################################
# Mouse mean
##################################################

mouse_mean <- df_long %>%
  filter(Group %in% c("AP","DP")) %>%
  group_by(
    Marker,
    Group,
    MouseID
  ) %>%
  summarise(
    Value = mean(Value, na.rm = TRUE),
    .groups = "drop"
  ) 
  

mouse_mean <- mouse_mean %>%
  filter(Group %in% c("AP","DP")) %>%
  mutate(
    Group = factor(
      Group,
      levels = c("AP","DP"),
      labels = c(
        "AM80 + αPD-L1",
        "DMSO + αPD-L1"
      )
    )
  )

##################################################
# Mixed model
##################################################
stats <- map_dfr(
  unique(df_long$Marker),
  function(marker){
    
    dat <- df_long %>%
      filter(
        Marker == marker,
        Group %in% c("AP","DP")
      )
    
    fit <- lmer(
      Value ~ Group +
        (1|MouseID),
      data = dat
    )
    
    emm <- emmeans(
      fit,
      ~ Group
    )
    
    p <- pairs(
      emm,
      adjust = "none"
    ) %>%
      as.data.frame() %>%
      pull(p.value)
    
    tibble(
      Marker = marker,
      p = p
    )
  }
)

##################################################
# Colors
##################################################
cols <- c(
  "AM80 + αPD-L1" = "#18499e",
  "DMSO + αPD-L1" = "#e7211a"
)
##################################################
# Theme
##################################################

theme_my <- theme_classic(base_size = 10) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(
      color = "black",
      size = 8
    )
  )

##################################################
# Plot function
##################################################
plot_marker <- function(marker_name){
  
  dat <- mouse_mean %>%
    filter(Marker == marker_name)
  
  pval <- stats %>%
    filter(Marker == marker_name) %>%
    pull(p)
  
  ymax <- max(dat$Value)
  
  ggplot(
    dat,
    aes(
      x = Group,
      y = Value
    )
  ) +
    
    geom_boxplot(
      aes(fill = Group),
      width = 0.30,
      alpha = 0.35,
      color = "black",
      outlier.shape = NA,
      linewidth = 0.8
    ) +
    
    geom_beeswarm(
      color = "black",
      size = 1,
      alpha = 0.4,
      cex = 5,
      priority = "density"
    ) +
    
    scale_fill_manual(
      values = cols
    ) +
    
    annotate(
      "segment",
      x = 1,
      xend = 2,
      y = ymax * 1.08,
      yend = ymax * 1.08
    ) +
    
    annotate(
      "text",
      x = 1.5,
      y = ymax * 1.13,
      label = paste0(
        "p = ",
        sprintf("%.3f", pval)
      ),
      size = 3
    ) +
    
    ggtitle(marker_name) +
    
    coord_cartesian(
      ylim = c(
        min(dat$Value),
        ymax * 1.2
      )
    ) +
    
    theme_my +
    theme(
      axis.text.x = element_text(
        angle = 30,
        hjust = 1
      )
    )
}

##################################################
# Create plots
##################################################

plist <- lapply(
  marker_info$Marker,
  plot_marker
)

##################################################
# Combine
##################################################

combined_plot <- wrap_plots(
  plist[c(5,6,7,8,9,10,4,3,1,2)],
  ncol = 3
)

combined_plot

##################################################
# Save
##################################################

ggsave(
  filename = "output/2026.06.13_AM80_PDL1_analysis.pdf",
  plot = combined_plot,
  width = 180,
  height = 180,
  units = "mm",
  device = cairo_pdf
)
