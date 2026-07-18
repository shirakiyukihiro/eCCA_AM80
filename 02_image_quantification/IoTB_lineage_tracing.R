# ============================================================
# Orthotopic (IoTB) model — lineage-tracing quantification (merge across images)
#
# Study      : Distinct cancer-associated fibroblast lineages and an
#              Am80-based stromal-reprogramming strategy in extrahepatic
#              cholangiocarcinoma (Meflin/ISLR).
# Repository : eCCA_AM80  (https://github.com/shirakiyukihiro/eCCA_AM80)
#
# Purpose    : Merge and summarise lineage-tracing quantifications in the orthotopic (implantation of tumour fragments; IoTB) model across imaged fields.
# Figure(s)  : Figure 4E (Col1a1, Pdgfrb, Acta2 among tdTomato+ Meflin-lineage cells; descriptive).
# Input      : data/IoTB_summary.xlsx.
# Output     : output/2026.06.10_IoTB_merge.pdf.
#
# Paths have been made relative for public release: place input files
# in ./data and write outputs to ./output (create these folders, or edit
# the paths below).
# ============================================================

library(readxl)
library(tidyverse)
library(patchwork)

file_path <- "data/IoTB_summary.xlsx"          

df <- read_excel(file_path, sheet = 1)

colnames(df)[1] <- "filename"
df$filename
df2 <- df%>% 
  mutate(
    marker = case_when(
      str_detect(filename, "ACTA2|Acta2") ~ "Acta2",
      str_detect(filename, "Col1a1") ~ "Col1a1",
      str_detect(filename, "Pdgfrb") ~ "Pdgfrb"
    ),
    percent = as.numeric(`%`),RFP = as.numeric(`RFP+`)
  )%>%mutate(ISH_RFP=merge/RFP*100)%>%dplyr::filter(!is.na(filename))

cols <- c(
  Acta2 = "#850000",
  Col1a1 = "#850000",
  Pdgfrb = "#850000"
)

plot_marker <- function(marker_name, fill_color){
  
  df_sub <- df2 %>%
    filter(marker == marker_name)
  
  summary_sub <- df_sub %>%
    summarise(
      mean = mean(ISH_RFP, na.rm = TRUE),
      sd   = sd(ISH_RFP, na.rm = TRUE)
    )
  
  ggplot(df_sub, aes(x = marker_name, y = ISH_RFP)) +
    
    geom_col(
      data = summary_sub,
      aes(x = marker_name, y = mean),
      width = 0.5,
      fill = fill_color,
      color = "black",
      linewidth = 0.6,
      inherit.aes = FALSE
    ) +
    
    geom_errorbar(
      data = summary_sub,
      aes(
        x = marker_name,
        ymin = mean - sd,
        ymax = mean + sd
      ),
      width = 0.12,
      linewidth = 0.7,
      inherit.aes = FALSE
    ) +
    
    geom_jitter(
      width = 0.08,
      height = 0,
      shape = 16,
      size = 2.5,
      color = "black",
      alpha = 0.8
    ) +
    
    theme_classic(base_size = 10) +
    
    labs(
      x = NULL,
      y = "Merge (%)"
    ) +scale_y_continuous(
      limits = c(0, 100)
    ) +
    
    theme(
      legend.position = "none",
      axis.title.x = element_blank()
    )
}

p_acta2 <- plot_marker("Acta2", cols["Acta2"])
p_col1a1 <- plot_marker("Col1a1", cols["Col1a1"])
p_pdgfrb <- plot_marker("Pdgfrb", cols["Pdgfrb"])

combined_plot <-
  p_acta2 +
  p_col1a1 +
  p_pdgfrb +
  plot_layout(ncol = 3)

combined_plot


ggsave(
  filename = "output/2026.06.10_IoTB_merge.pdf",
  plot = combined_plot,
  width = 90,
  height = 50,
  units = "mm",
  device = cairo_pdf
)
