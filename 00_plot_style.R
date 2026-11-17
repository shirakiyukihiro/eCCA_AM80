## =====================================================================
## 00_plot_style.R  —  shared house style for all quantification figures
## Source this at the top of each figure script:  source("00_plot_style.R")
## Style basis = Fig 6/7/8 (box + beeswarm dots + p annotation),
## with the y-axis title shown so script output matches the final figures.
##
## v3: balanced look - both box and dots semi-transparent. Box keeps a
##     medium-grey outline (grey35) with a light fill so it is not buried;
##     dots are pch 21 with a THIN white edge and ~0.55-0.6 alpha so they
##     read clearly without overpowering the box.
## =====================================================================

library(ggplot2)
library(ggbeeswarm)

## fixed shape palette (optional shape-by-individual mode; up to 10)
shape_pal <- c(16, 17, 15, 18, 3, 4, 8, 7, 5, 6)

## colour-blind-safe palette for per-individual dots (Okabe-Ito; up to 8)
mouse_pal <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7", "#999999")

## ---- unified theme -------------------------------------------------
theme_my <- theme_classic(base_size = 10) +
  theme(
    legend.position = "none",
    axis.title.x    = element_blank(),
    axis.title.y    = element_text(size = 9, colour = "black"),
    axis.text       = element_text(size = 9, colour = "black"),
    axis.line       = element_line(colour = "black"),
    plot.title      = element_text(hjust = 0.5, size = 10)
  )

## ---- unified plotting helper --------------------------------------
## dat     : data frame whose ROWS are the points drawn as dots
##           (per-mouse means for per-animal graphs; individual fields
##            for image-unit graphs)
## xvar/yvar/fillvar : column names (strings)
## ylab    : y-axis title (string or expression())
## pval    : numeric p to annotate, or NULL (descriptive, no p)
## cols    : named colour vector for scale_fill_manual, or NULL
## xlabels : named vector to relabel x ticks, or ggplot2::waiver()
## p_x     : x position of the p label (default 1.5 = between 2 groups)
## ylim    : fixed y-axis range, e.g. c(0, 100); NULL = auto-scale
## idvar   : column name to colour dots by individual (per-animal mode);
##           NULL = single haloed black dots (default)
plot_box <- function(dat, xvar, yvar, fillvar, ylab,
                     pval = NULL, cols = NULL,
                     xlabels = ggplot2::waiver(), p_x = 1.5, ylim = NULL,
                     idvar = NULL) {

  ymax <- max(dat[[yvar]], na.rm = TRUE)
  ymin <- min(dat[[yvar]], na.rm = TRUE)

  g <- ggplot(dat, aes(x = .data[[xvar]], y = .data[[yvar]]))

  if (is.null(idvar)) {
    ## ---- default: marker-coloured box + white-haloed black dots ----
    g <- g +
      geom_boxplot(
        aes(fill = .data[[fillvar]]),
        width = 0.50, alpha = 0.35, colour = "grey35",
        linewidth = 0.4, outlier.shape = NA, staplewidth = 0.5
      ) +
      geom_beeswarm(
        shape = 21, fill = "black", colour = "white",
        stroke = 0.25, size = 1.7, alpha = 0.55, cex = 3
      )
    if (!is.null(cols)) g <- g + scale_fill_manual(values = cols)

  } else {
    ## ---- per-individual: neutral open box + colour-by-mouse dots ---
    ## box drawn unfilled with a thin grey outline so nothing competes
    ## with the dots; each dot is a white-bordered, colour-filled circle
    g <- g +
      geom_boxplot(
        fill = "grey50", alpha = 0.18, colour = "grey35",
        width = 0.50, linewidth = 0.4, outlier.shape = NA, staplewidth = 0.5
      ) +
      geom_beeswarm(
        aes(fill = .data[[idvar]]),
        shape = 21, colour = "white", stroke = 0.25,
        size = 1.8, alpha = 0.6, cex = 3
      ) +
      scale_fill_manual(values = mouse_pal, name = "Mouse", drop = FALSE)
  }

  g <- g +
    scale_x_discrete(labels = xlabels) +
    labs(x = NULL, y = ylab) +
    theme_my

  ## per-individual graphs need the legend that theme_my switches off
  if (!is.null(idvar)) {
    g <- g + theme(
      legend.position = "bottom",
      legend.title    = element_text(size = 8, colour = "black"),
      legend.text     = element_text(size = 8, colour = "black"),
      legend.key.size = unit(3, "mm")
    )
  }

  ## y-range + p-label position
  if (!is.null(ylim)) {
    g   <- g + coord_cartesian(ylim = ylim)
    p_y <- ylim[2] * 0.96
  } else {
    g   <- g + coord_cartesian(ylim = c(min(0, ymin), ymax * 1.18))
    p_y <- ymax * 1.08
  }

  if (!is.null(pval)) {
    lab <- if (pval < 0.001) "p < 0.001" else paste0("p = ", signif(pval, 2))
    g <- g + annotate("text", x = p_x, y = p_y, label = lab, size = 3)
  }
  g
}

## ---- shared y-axis labels (marker -> axis title) ------------------
ylab_map <- c(
  MT      = "Fibrosis area (%)",
  Pdgfrb  = "PDGFR\u03b2+ area (%)",
  PDGFRB  = "PDGFR\u03b2+ area (%)",
  SMA     = "\u03b1-SMA+ area (%)",
  Ki67    = "Ki-67+ cells (%)",
  CD31    = "CD31+ vascular area (%)",
  CD11b   = "CD11b+ area (%)",
  CD86    = "CD86+ area (%)",
  CD3     = "CD3+ cells (/HPF)",
  CD4     = "CD4+ cells (/HPF)",
  CD8     = "CD8+ cells (/HPF)",
  CD163   = "CD163+ cells (/HPF)"
)
