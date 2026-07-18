# ============================================================
# Subcutaneous tumour IHC quantification — Am80 + chemotherapy
#
# Study      : Distinct cancer-associated fibroblast lineages and an
#              Am80-based stromal-reprogramming strategy in extrahepatic
#              cholangiocarcinoma (Meflin/ISLR).
# Repository : eCCA_AM80  (https://github.com/shirakiyukihiro/eCCA_AM80)
#
# Purpose    : Read subcutaneous-tumour immunostaining quantifications for the Am80 + chemotherapy arms (LMM-based analysis and plots). Shares the input workbook with the WT-vs-KO and Am80/PD-L1 IHC scripts.
# Figure(s)  : Figure 7E-J (subcutaneous IHC box plots, Am80 + chemotherapy).
# Input      : data/subcutaneous_IHC_data.xlsx.
# Output     : output/2026.06.13_AM80Chemo_subcutaneous_IHC.pdf.
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
  sheet = sheets[2],
  skip = 1
)

marker_info <- tibble(
  
  Marker = c(
    "Ki67",
    "CD31",
    "SMA",
    "MT",
    "PDGFRB",
    "CD3",
    "CD11b"
  ),
  
  WTAC = c(
    "WT-AC...1",
    "WT-AC...8",
    "WT-AC...15",
    "WT-AC...22",
    "WT-AC...29",
    "WT-AC...36",
    "WT-AC...44"
  ),
  
  WTDC = c(
    "WT-DC...2",
    "WT-DC...9",
    "WT-DC...16",
    "WT-DC...23",
    "WT-DC...30",
    "WT-DC...37",
    "WT-DC...45"
  ),
  
  KOAC = c(
    "KO-AC...3",
    "KO-AC...10",
    "KO-AC...17",
    "KO-AC...24",
    "KO-AC...31",
    "KO-AC...38",
    "KO-AC...46"
  ),
  
  KODC = c(
    "KO-DC...4",
    "KO-DC...11",
    "KO-DC...18",
    "KO-DC...25",
    "KO-DC...32",
    "KO-DC...39",
    "KO-DC...47"
  )
)

group_structure <- c(
  "WT-AC" = 5,
  "WT-DC" = 5,
  "KO-AC" = 6,
  "KO-DC" = 4
)

make_group_df <- function(values,
                          marker,
                          group){
  
  values <- values[!is.na(values)]
  
  n_mouse <- group_structure[group]
  
  mouse_id <- rep(
    paste0(group, "_", seq_len(n_mouse)),
    each = 5
  )
  
  tibble(
    Marker = marker,
    Group = group,
    MouseID = mouse_id[seq_along(values)],
    Value = values
  )
}

df_long <- map_dfr(
  seq_len(nrow(marker_info)),
  function(i){
    
    marker <- marker_info$Marker[i]
    
    bind_rows(
      
      make_group_df(
        df[[marker_info$WTAC[i]]],
        marker,
        "WT-AC"
      ),
      
      make_group_df(
        df[[marker_info$WTDC[i]]],
        marker,
        "WT-DC"
      ),
      
      make_group_df(
        df[[marker_info$KOAC[i]]],
        marker,
        "KO-AC"
      ),
      
      make_group_df(
        df[[marker_info$KODC[i]]],
        marker,
        "KO-DC"
      )
      
    )
    
  }
)

df_long_model <- df_long %>%
  separate(
    Group,
    into = c("Genotype","Region"),
    sep = "-"
  )

mouse_mean <- df_long %>%
  group_by(
    Marker,
    Group,
    MouseID
  ) %>%
  summarise(
    Value = mean(Value),
    .groups = "drop"
  ) %>%
  separate(
    Group,
    into = c("Genotype","Region"),
    sep = "-"
  )

stats <- map_dfr(
  unique(df_long_model$Marker),
  function(marker){
    
    dat <- df_long_model %>%
      filter(Marker == marker)
    
    fit <- lmer(
      Value ~ Genotype * Region +
        (1|MouseID),
      data = dat
    )
    
    emm <- emmeans(
      fit,
      ~ Region | Genotype
    )
    
    pairs(
      emm,
      adjust = "holm"
    ) %>%
      as.data.frame() %>%
      mutate(
        Marker = marker
      )
    
  }
)

stats

stats <- stats %>%
  mutate(
    label = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )


mouse_mean <- mouse_mean %>%
  mutate(
    Group = factor(
      paste(
        Genotype,
        Region,
        sep = "-"
      ),
      levels = c(
        "WT-AC",
        "WT-DC",
        "KO-AC",
        "KO-DC"
      )
    )
  )

cols <- c(
  "WT-AC" = "#18499e",
  "WT-DC" = "#e7211a",
  "KO-AC" = "#6e8fd4",
  "KO-DC" = "#f08c86"
)


theme_my <- theme_classic(base_size = 10) +
  theme(
    legend.position = "none",
    axis.title = element_blank(),
    axis.text = element_text(
      color = "black",
      size = 8
    )
  )

plot_marker <- function(marker_name){
  
  dat <- mouse_mean %>%
    filter(Marker == marker_name)
  
  pdat <- stats %>%
    filter(Marker == marker_name)
  
  ymax <- max(dat$Value)
  
  wt_lab <- pdat %>%
    filter(Genotype == "WT") %>%
    pull(p.value)%>%round( 3)
  
  ko_lab <- pdat %>%
    filter(Genotype == "KO") %>%
    pull(p.value)%>%round( 3)
  
  ggplot(
    dat,
    aes(
      x = Group,
      y = Value
    )
  ) +
    
    geom_boxplot(
      aes(fill = Group),
      width = 0.5,
      alpha = 0.6,
      color = "black",
      outlier.shape = NA
    ) +
    
    geom_beeswarm(
      color = "black",
      size = 1.3,
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
      y = ymax*1.05,
      yend = ymax*1.05
    ) +
    
    annotate(
      "text",
      x = 1.5,
      y = ymax*1.10,
      label = wt_lab,
      size = 4
    ) +
    
    annotate(
      "segment",
      x = 3,
      xend = 4,
      y = ymax*1.05,
      yend = ymax*1.05
    ) +
    
    annotate(
      "text",
      x = 3.5,
      y = ymax*1.10,
      label = ko_lab,
      size = 4
    ) +
    
    coord_cartesian(
      ylim = c(
        min(dat$Value),
        ymax*1.20
      )
    ) +
    
    ggtitle(marker_name) +
    
    theme_my
}

p1 <- plot_marker("Ki67")
p2 <- plot_marker("CD31")
p3 <- plot_marker("SMA")
p4 <- plot_marker("MT")
p5 <- plot_marker("PDGFRB")
p6 <- plot_marker("CD3")
p7 <- plot_marker("CD11b")

combined_plot <-
  (p4 + p3 + p1) /
  (p2 + p6 + p7) 

combined_plot

ggsave(
  filename = "output/2026.06.13_AM80Chemo_subcutaneous_IHC.pdf",
  plot = combined_plot,
  width = 180,
  height = 120,
  units = "mm",
  device = cairo_pdf
)

