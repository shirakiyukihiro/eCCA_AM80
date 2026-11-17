# =============================================================================
# Figure panels for the mouse biliary fibroblast single-cell analysis
#
# Study    : Identification of distinct Meflin+ fibroblast subsets in the
#            biliary stroma reveals a therapeutic strategy for extrahepatic
#            cholangiocarcinoma
# Figures  : Figure 3, Supplementary Figure S1
#
# This is the submission copy of the panel code. It redraws the panels from
# the objects written by 01_mouse_biliary_fibroblast_scRNAseq.R and writes to
# CONFIG$out_dir/figures_final/. No analysis is repeated here.
#
# FIGURE 3, published panel order
#   A1  whole-tissue UMAP, fibroblasts highlighted
#   A2  fibroblast UMAP, three subsets
#   B   marker dot plot
#   C   feature plots, 6 genes
#   D   Gene Ontology, one versus rest
#   E   GSEA normalised enrichment scores, one versus rest
#   F   module scores, universal fibroblast programmes
#
# The order is the argument: clustering gives three subsets (A), their markers
# define them (B), a subset of those markers locates them on the UMAP (C), GO
# characterises what each subset is doing (D), and the module scores place the
# subsets against the published universal-fibroblast types (F), which is where
# the paper hands over to the lineage tracing in Figures 4 and 5.
#
# Redraws every figure from the saved objects. No analysis is repeated, so
# adjusting a font size or a point size costs seconds instead of an hour.
#
#   source("02_figure_panels_Fig3_S1.R")  # CONFIG, helpers, theme, functions
#   fig_all()                  # everything
#   fig_mouse()                # the published panels
#   fig_umaps(); fig_dotplots(); fig_violins(); fig_go()   # one kind
#
# Everything the figures need comes from disk:
#   02_allcells_annotated.rds        all mouse cells, annotated
#   03c_fibroblast_subsets.rds       mouse fibroblasts with subsets
#   05_go_results.rds                enrichGO results (section 13)
#
# The human-dataset branches (fig_human() and the explor_* helpers) are inert
# unless CONFIG$human_rds is set; that analysis is not part of the published
# study and the objects are not distributed.
#
# Run the analysis script at least once to produce them.
# =============================================================================

# Printed on load, together with the absolute path the file was read from. The
# point is to make "I edited the file but nothing changed" impossible to miss:
# if this string is not the one in the message that came with the file, R is
# sourcing a different copy.
FIG_VERSION <- "final-2026-10-08 (Fig 3E GSEA NES; Fig 3F module scores; Fig S1D named genes)"

# CONFIG and the helper functions come from sections 0 and 1 of the analysis
# script, so there is one definition of the paths and of avg_expr(), fig(),
# obj() and the rest. Nothing else from that script is executed.
SCRIPT <- "mouse_biliary_fibroblast_scRNAseq.R"
local({
  src <- readLines(SCRIPT, warn = FALSE)
  at  <- grep("^#@ STEP [0-9]+$", src)
  num <- as.integer(sub("^#@ STEP ", "", src[at]))
  end <- c(at[-1] - 1L, length(src))
  for (s in c(0L, 1L)) {
    k <- which(num == s)
    eval(parse(text = paste(src[at[k]:end[k]], collapse = "\n")),
         envir = globalenv())
  }
})


# =============================================================================
# FIGURE SETTINGS
# -----------------------------------------------------------------------------
# All the knobs in one place. Change a value here and re-run one fig_* function.
# =============================================================================

FIG <- list(
  # Subdirectory of CONFIG$out_dir that save_plot() writes into. "figures" is
  # what figures.R and the analysis script's fig() use; this script uses its
  # own directory so that re-running either of those cannot overwrite a panel
  # that has already been placed in the manuscript figure. Do not point this
  # back at "figures".
  out_subdir = "figures_final",

  # Journals ask for Arial or Helvetica. On Linux "Arial" exists only if the
  # msttcorefonts are installed; setup_font() below reports what was actually
  # found rather than letting the device substitute silently, which is how a
  # figure ends up in DejaVu Sans without anyone noticing.
  font_preference = c("Arial", "Helvetica", "Nimbus Sans", "Liberation Sans",
                      "sans"),

  # CHANGED 2026-09-30: 7 -> 8 pt for axis and legend text, 8 -> 9 for panel
  # titles. Everything downstream is expressed as base_size or base_size + 1,
  # so this is the only place to change it -- but panels whose PANEL_SIZE was
  # tuned at 7 pt now hold about 14 % more text. PANEL_SIZE$go is retuned below;
  # if another panel comes out crowded, the width, not the font, is what to
  # change back.
  # ADDED 2026-09-30. Every panel is written twice: a PDF (as before) and an
  # SVG. Open the SVG in Illustrator, not the PDF. Two reasons, both of which
  # were costing time on every figure:
  #   1. cairo_pdf positions text run by run, so one axis label arrives in
  #      Illustrator as several separate point-text objects that come apart
  #      when the artwork is resized. svglite writes one <text> element per
  #      label, so each label is one editable object.
  #   2. A PDF embeds the glyphs of whatever font the system actually had, so
  #      the file says "Nimbus Sans" and has to be converted to Arial by hand.
  #      An SVG only names the family, so Illustrator resolves "Arial" against
  #      the real Arial installed on the machine that opens it.
  # Needs the svglite package; without it only the PDF is written and a message
  # says so.
  also_svg    = TRUE,

  base_size   = 8,      # point size for axis text; titles scale from it
  title_size  = 9,
  strip_size  = 8,

  # Rasterisation. Points become a raster; axes, text and lines stay vector, so
  # the file opens quickly in Illustrator and the labels remain editable.
  raster     = TRUE,
  raster_px  = c(2000, 2000),   # SEE BELOW -- this is not a dpi
  raster_dpi = 600,             # used only by ggrastr, where it IS a dpi

  # WHY THE POINTS DISAPPEARED
  #
  # With raster = TRUE, Seurat does not draw ggplot points at all: it hands the
  # data to scattermore, where the canvas size is given in PIXELS and pt.size is
  # the point radius in PIXELS, not in millimetres. The previous settings drew
  # one-pixel points on a 600 x 600 canvas; stretched to a 4.6-inch figure that
  # is 130 dpi with dots too small to see, and in a feature plot the coloured
  # cells vanish under the grey ones because each is a single pixel.
  #
  # Both numbers therefore have to move together. The canvas must be large
  # enough for the printed size (2000 px across 4.6 in is 430 dpi) and the point
  # radius must be in the same units: 4 px on that canvas is about 0.7 pt, which
  # is a visible dot.
  #
  # To tune: change umap_point first. If the plot turns into a solid mass,
  # lower it; if it is still faint, raise it. Only raise raster_px if the dots
  # start to look square. raster = FALSE draws true vector points and is the way
  # to check what the raster is costing you, at the price of a heavy file.
  # Point size is given for a REFERENCE number of cells and scaled from there by
  # pt_for() below, because one value cannot suit every panel: the mouse
  # whole-tissue plot has 3,700 cells, the fibroblast subset 1,800 and the human
  # CAF subset 364. A size that reads well at 3,700 is invisible at 364. The
  # scaling is by sqrt(reference / n), which keeps the fraction of the canvas
  # covered by ink roughly constant.
  umap_point       = 8,      # radius in px at umap_point_ref cells
  umap_point_ref   = 3700,
  # ...AND at a panel this wide, in inches. scattermore draws onto a fixed
  # 2000 x 2000 pixel canvas and pt.size is a radius in canvas pixels, so the
  # printed dot diameter is proportional to the panel width: the same value
  # gives dots twice as large on a panel twice as wide. Every plotting function
  # takes a width_in argument and scales the radius by width_in / this.
  #
  # This replaces the composite-specific overrides that used to exist. Those
  # were a second set of point sizes that had to be kept in step with the first
  # by hand, and they were not; a value tuned on a standalone panel then looked
  # wrong in the composite and vice versa.
  # THESE TWO NUMBERS ARE THE ONLY LEVER YOU SHOULD TOUCH for dot size.
  # Each is the panel width, in inches, at which the point size above was judged
  # to look right. Everything else follows: pt_for() scales the radius by
  # (actual width / this). Set them wrong and every panel is wrong by the same
  # factor -- which is what happened when feature plots were scaled against the
  # UMAP reference: feature_point = 8 was tuned on a 2.2-inch SUB-PLOT, not on a
  # 4.6-inch panel, so Figure 5C came out at a third of the intended size and
  # the expressing cells disappeared.
  #
  # umap_point_ref_in  : the width of a whole UMAP panel.
  # feature_point_ref_in: the width of ONE GENE'S sub-plot, not of the row.
  umap_point_ref_in    = 6.0,
  feature_point_ref_in = 2.2,
  umap_point_range = c(3, 20),   # never smaller or larger than this
  umap_alpha       = 0.9,

  # Superscripts in the cluster labels.
  #
  # The labels are stored plainly -- "Meflin+/Thy1-high" -- so that they stay
  # legible in the object and in every table. Here they are turned into
  # "Meflin<sup>+</sup>/Thy1<sup>high</sup>" and drawn by ggtext, which renders
  # real superscripts as vector text.
  #
  # Unicode was the alternative and was rejected: there is no superscript "high"
  # that a journal-safe font carries. U+1D4D and friends exist but Arial,
  # Helvetica and the Nimbus clones do not have them, so the characters would
  # silently become boxes in the PDF. ggtext has a dependency; missing glyphs
  # have no fix at all.
  #
  # Set FALSE to draw the labels exactly as stored.
  rich_labels      = TRUE,

  # Feature plots get their own floor. The scaling above is right for a DimPlot,
  # where every point carries a colour and a dense panel is still readable, but
  # wrong for a feature plot of a rare population: in a 30,000-cell whole-tissue
  # object the fibroblasts are a few hundred cells, and at the 2 px floor they
  # are the only thing worth seeing and the only thing invisible. The cost of a
  # larger floor is that the grey majority gets denser, which does not matter
  # because the grey majority is not the message.
  feature_point       = 9,
  feature_point_range = c(4, 20),

  # Feature plot colour scale. Expression is zero in most cells and high in a
  # few, so a linear scale over the full range paints almost everything the low
  # colour and the signal disappears.
  #   "expressing"  per gene, 0 to the given quantile of the NON-ZERO values
  #   "all"         Seurat's own q5/q95 over every cell
  #   NA            the full range
  # "expressing" is the default because a quantile taken over every cell breaks
  # exactly where it matters: if a gene is detected in 3 % of cells, the 95th
  # percentile over all cells is 0, the scale collapses and the panel comes out
  # blank. Taking the quantile among the cells that express the gene puts the
  # contrast inside the positive population, whatever fraction of the tissue it
  # is. The top few per cent saturate, which is intended.
  feature_scale        = "expressing",
  # Lower than 0.95: with the top 5 % setting the saturation point, the bulk of
  # expressing cells were compressed into the pale half of the ramp.
  feature_hi_quantile  = 0.90,
  feature_min_express  = 10,   # fewer positive cells than this: use their max
  feature_min_cutoff   = "q5",  # used only when feature_scale = "all"
  feature_max_cutoff   = "q95",
  # Two stops keep zero grey and the maximum red, but a weakly expressing cell
  # then sits next to grey and is hard to pick out. Adding a middle stop, e.g.
  #   c("grey90", "#F0A868", "#C1443C")
  # makes any detected cell visibly warm at the price of a less linear ramp.
  # Four stops, not two. With a straight grey-to-red ramp a cell at a fifth of
  # the range is a barely tinted grey, and most detected cells sit in the lower
  # part of the range -- so "detected" and "not detected" looked the same. The
  # stops are evenly spaced by scale_color_gradientn, so putting most of the
  # colour change in the first third means any expressing cell reads as red.
  # Zero still maps to grey.
  feature_cols       = c("grey92", "#E79B92", "#D4685C", "#C1443C"),

  # Cluster labels drawn on top of the feature plots. Without them the reader has
  # to hold the UMAP layout in memory while looking at six panels.
  #   "auto"    names when there are few short ones, numbers otherwise
  #   "number"  a numeral per cluster, with the key in the companion DimPlot
  #   "name"    the cluster label itself
  #   "none"    no labels
  # "auto" is the default because the right answer differs by panel: three
  # clusters called Thy1-high and Cxcl12-high fit on a small panel,
  # twelve cell-type names repeated six times do not.
  # "number", not "auto". auto compares the subset names against
  # feature_label_max_ch = 20; they are 17 characters, so it was choosing names
  # -- three twenty-character labels with superscripts, laid over the cells they
  # point at, on a panel a couple of inches wide. The numerals are keyed by the
  # UMAP panel beside them and are drawn with a halo, so they cover nothing.
  # Set "name" per call if a standalone panel needs to stand on its own.
  feature_label        = "number",
  feature_label_size   = 2.0,

  # How a cluster label is kept legible over the points underneath.
  #
  #   "halo"  the numeral is drawn eight times in white, offset by a fraction of
  #           a millimetre, with the black glyph on top. Nothing is covered
  #           except the few pixels the strokes themselves occupy.
  #   "box"   a translucent white rectangle behind the text. On a feature plot
  #           this sits exactly over the cells that express the gene -- the
  #           hottest, reddest part of the panel, which is the one place a label
  #           must not hide. That is what it was doing.
  #
  # Name labels carry superscripts and need geom_richtext, which cannot be
  # haloed, so they keep the box. Numerals have no markup and use the halo.
  label_style          = "halo",
  label_halo_r         = 0.12,   # shadowtext outline radius
  label_halo_frac      = 0.012,  # fallback offset, fraction of the embedding
  feature_label_max_n  = 4,    # at most this many clusters for names
  feature_label_max_ch = 20,   # and labels no longer than this

  # How the companion DimPlot carries the key for numbered feature plots.
  #   "name"    names repelled onto the plot, no legend (the original)
  #   "number"  the same numerals as the feature plots on the plot, and a
  #             legend reading "1  Fibroblast-1", "2  Fibroblast-2", ...
  #   "auto"    whichever the companion feature plots resolved to
  # "number" removes the need to transcribe a key into the figure caption. The
  # numbers are positions in the plotting order, not Seurat cluster ids: the
  # legend therefore always reads 1, 2, 3, ... downwards. Seurat's own cluster
  # numbers are not used because they are an artefact of the run (they move with
  # resolution and seed) and are not reported anywhere in the manuscript.
  # "auto" is the default so that the three-subset panels, whose feature plots
  # are labelled with names, keep names on the DimPlot as well: a numbered key
  # is only worth its legend when the names are too long to sit on the panels.
  umap_key             = "auto",
  umap_key_gap         = "  ",  # between numeral and name in the legend
  umap_label_size      = 2.6,   # numerals on the DimPlot itself

  # Jittered points over violins and boxes. The problem with overlaying points
  # on a thousand cells is that the layer goes solid black and hides the
  # distribution it was meant to show. Three controls fix that, and they have to
  # be used together: draw a bounded number of points per group, make them small,
  # and make them transparent.
  # Start here and adjust in this order: size until single points are legible,
  # then max_n until the cloud stops filling the violin, then alpha for contrast.
  # Raised from 0.6 / 0.35: at eleven facets the points were being drawn small
  # AND faint, and the two compound into invisible. Size first, then alpha; drop
  # jitter_max_n before lowering either if the cloud starts filling the violin.
  jitter_max_n     = 350,    # points drawn per group and panel; NA draws all
  jitter_point     = 0.6,
  jitter_alpha     = 0.55,
  jitter_width     = 0.22,
  jitter_colour    = "grey15",

  # Box plots: the box is drawn semi-transparent so points underneath show
  # through, and outliers are suppressed because the jitter layer already shows
  # every point.
  box_alpha        = 0.55,
  box_linewidth    = 0.3,
  violin_linewidth = 0.3,

  # Marker dot plots. 45 degrees is the convention for gene names on an x axis:
  # vertical labels cost height and are harder to read at small point sizes.
  dot_angle        = 45,

  # Enrichment dot plots
  # Enrichment dot plots. Three channels, and the point is that each should
  # carry something different.
  #   go_x      "fold"   fold enrichment, (k/n) / (M/N): the effect size
  #             "ratio"  gene ratio, k/n
  #   go_order  "x"      terms ordered by whatever is on the x axis
  #             "p"      terms ordered by adjusted P
  # Defaults are "fold" and "x". Gene ratio is k/n with n fixed within a panel,
  # so it is proportional to Count and duplicates the dot size; fold enrichment
  # is not. Ordering by adjusted P while plotting something else on x makes the
  # dots zigzag and also spends the y axis and the colour on the same variable.
  # Terms are still SELECTED by adjusted P -- only the display order changes --
  # which has to be stated in the legend.
  # Fixed aspect for UMAP-derived panels. 1 keeps the embedding undistorted and
  # makes every panel drawn on it the same shape. NULL lets the panel stretch.
  umap_aspect      = 1,

  # Numerals on the plot with the names in the key beside it: the subset names
  # carry superscripts and run to about twenty characters, which covers the
  # cells they are meant to label.
  fib_umap_key     = "number",

  # Which whole-tissue clusters the highlighted UMAP colours, and in what.
  # One colour, not three: the whole-tissue clusters Fibroblast-1/2/3 do NOT
  # correspond one-to-one to the three subsets (145 of the 460 Fibroblast-2
  # cells end up in the Cxcl12-high subset), so three colours beside three
  # subsets would invite the reader to pair them. One colour claims nothing.
  highlight_keep   = c("Fibroblast-1", "Fibroblast-2", "Fibroblast-3"),
  highlight_grey   = "grey85",
  highlight_colour = NULL,          # NULL = palette[1]
  # "number" numbers every whole-tissue cluster, keyed by Fig. S1A;
  # "keep" numbers only the extracted ones; "none" draws no labels.
  highlight_label  = "number",
  # Printed dot diameter for the highlighted panel, as a fraction of the
  # feature-plot dot. Not 1: this panel holds 5,105 cells in 2.6 in against
  # 1,784 in 2.4 in, so equal dots would be a solid mass. At the old setting it
  # was 0.23 of the feature-plot dot and had hit the lower clamp, which is why
  # nothing on the panel was visible.
  highlight_dot_frac = 0.6,

  # CHANGED 2026-09-29: Col1a1 and Pi16 added, 4 genes -> 6, laid out 3 x 2.
  # The rows are the panel's argument, so keep this order:
  #   row 1  Islr, Col1a1, Pi16  -- shared. Islr and Col1a1 say all three
  #          clusters are fibroblasts; Pi16 is present in two of the three and
  #          absent from the Cxcl12-high subset, which is the pattern the human
  #          data repeat, and it is here as a finding, not as a name.
  #   row 2  Thy1, Mfap5, Cxcl12 -- one discriminator per subset, in subset
  #          order, so the row reads left to right as PBF, PPF, CxF.
  fib_feature_genes    = c("Islr", "Col1a1", "Pi16",
                           "Thy1", "Mfap5", "Cxcl12"),
  # 3 x 2 rather than 6 in a row or 2 x 3. The sub-plot width is what sets the
  # printed dot radius (see fig_panels), and 7.2 / 3 keeps it at the 2.4 in the
  # 2 x 2 layout had, so Fig3A2 stays matched without retuning. A 2 x 3 would
  # be 7.9 in tall, which is taller than the page allows beside Fig3A.
  fib_feature_ncol     = 3,
  fib_feature_genes_s1 = c("Islr", "Thy1", "Pi16", "Cxcl12", "Col1a1", "Acta2"),
  allcell_feature_genes = c("Islr", "Thy1", "Col1a1", "Pdgfra", "Pdgfrb",
                            "Acta2"),

  # NULL leaves whatever Seurat derived from the reduction key.
  umap_axis_titles = c("UMAP1", "UMAP2"),

  go_x             = "fold",
  go_order         = "x",
  # Applied before the top-n selection. A term of three genes can reach a fold
  # enrichment of 15 and would otherwise head the plot on the strength of three
  # cells' worth of signal.
  go_min_count     = 10,
  # Maximum size of the term in the background annotation. GO-BP contains very
  # broad terms -- "ossification", "heart morphogenesis", "visual system
  # development" -- which are annotated to hundreds of genes, overlap heavily
  # with any ECM or developmental programme, and reach significance in almost
  # any fibroblast contrast. They are not wrong, they are uninformative. NA
  # keeps every term; 500 is a common cut, 300 is aggressive.
  #
  # Set back to NA after testing: 500 and even 300 removed nothing from either
  # panel. The generic terms in the Thy1-high panel -- osteoblast
  # differentiation (183 genes), heart morphogenesis (225), reproductive
  # structure development (265) -- are not large, they are simply unrelated to
  # this tissue, and term size does not measure that. A filter that changes no
  # result still has to be described in the Methods, so it is off.
  go_max_bg        = NA,
  # Terms always kept, whatever the ranking, given as exact GO descriptions.
  # Use only for terms named as hypotheses BEFORE the analysis, and say so in
  # the legend; otherwise this is simply selecting the result you wanted.
  go_include       = NULL,
  go_n_terms       = 12,
  # CHANGED 2026-09-29: 24 -> 28 -> 32. Wrapping later puts fewer lines under
  # each term, which is what stops the names touching. The second increase pays
  # for the height being cut to three quarters: at 32 characters most GO-BP
  # names fit on two lines instead of three, so the same twelve terms fit in
  # 4.5 in. It also widens the term column, but the panel got 10 % wider at the
  # same time, so the plotting area does not lose anything.
  #
  # CHANGED 2026-09-30: 32 -> 24. Two things pull the panel wider at once --
  # the font went from 7 to 8 pt, and the width has to come down. Both are paid
  # for out of the term column, which is the widest part of the panel: at 24
  # characters and 8 pt a term name is about 1.3 in instead of 1.6 in, three
  # times over. The cost is more three-line names, which is why the height goes
  # up rather than down.
  go_wrap          = 24,     # characters before a term name wraps
  # CHANGED 2026-09-29: "column" -> "row". Stacked, the panel had to fit
  # 3 x 12 = 36 wrapped term names into one column of 6.2 in, about 0.17 in per
  # term, and multi-line names overlapped. Side by side, each column carries
  # 12 terms over the full height, so the same terms get roughly three times
  # the vertical space without the panel growing at all in that direction.
  go_layout        = "row",      # the three subsets side by side, shared guides
  go_legend        = "bottom",   # "right" only pays off on a wide panel
  go_dot_range     = c(1, 4),

  palette = c("#C1443C", "#3C6EC1", "#4E9A51", "#8A5FA8", "#D98C2B",
              "#5FA8A0", "#B4656F", "#7D7D7D")
)

suppressPackageStartupMessages({
  library(Seurat); library(tidyverse); library(patchwork)
})
has_pkg <- function(p) requireNamespace(p, quietly = TRUE)
if (!has_pkg("ggrastr"))
  message("ggrastr is not installed; violin and box point layers will be vector. ",
          "install.packages('ggrastr') to rasterise them.")

# --- font ---------------------------------------------------------------------
setup_font <- function(prefs = FIG$font_preference) {
  available <- tryCatch({
    if (has_pkg("systemfonts")) unique(systemfonts::system_fonts()$family)
    else character(0)
  }, error = function(e) character(0))

  if (!length(available)) {
    message("Cannot enumerate system fonts (install 'systemfonts' to check). ",
            "Using '", prefs[1], "' and trusting the device to find it.")
    return(prefs[1])
  }
  hit <- prefs[prefs %in% available]
  if (length(hit) && identical(hit[1], prefs[1])) return(prefs[1])

  # CHANGED 2026-09-30. The preferred family is not installed under its own
  # name. Before falling back, ask fontconfig what it would substitute: on
  # essentially every Linux distribution "Arial" is aliased to Liberation Sans
  # or Nimbus Sans, which are metric clones, so R already draws Arial widths --
  # it simply records the clone's name in the file.
  #
  # When that is what is happening, keep the requested NAME. The glyphs are
  # unchanged (fontconfig returns the same file either way) and every
  # measurement is unchanged, but the SVG now carries font-family="Arial", so
  # Illustrator opens it in the real Arial instead of asking for a conversion
  # on every file. This is only safe for a metric clone: naming a font whose
  # widths differ from the ones R measured would move every line break, so an
  # unrecognised substitute falls through to the old behaviour.
  metric_clones <- c("Liberation Sans", "Nimbus Sans", "Nimbus Sans L",
                     "Arimo", "Helvetica", "Albany AMT", "FreeSans")
  subst <- tryCatch({
    if (!has_pkg("systemfonts")) NA_character_ else {
      pth <- systemfonts::match_font(prefs[1])$path
      sf  <- systemfonts::system_fonts()
      sf$family[match(pth, sf$path)]
    }
  }, error = function(e) NA_character_)
  if (!is.na(subst) && subst %in% metric_clones) {
    message("'", prefs[1], "' is not installed, but fontconfig substitutes '",
            subst, "', which is metrically compatible. Keeping the name '",
            prefs[1], "' so the SVG opens in Arial without conversion.")
    return(prefs[1])
  }

  if (!length(hit)) {
    warning("None of ", paste(prefs, collapse = ", "), " is installed and the ",
            "substitute for '", prefs[1], "' is ",
            if (is.na(subst)) "unknown" else paste0("'", subst, "'"),
            ", which is not a known metric clone. The figures will not be in ",
            "the family you asked for. Falling back to 'sans'.")
    return("sans")
  }
  message("'", prefs[1], "' not installed; using '", hit[1], "'. ",
          "Helvetica and the Nimbus/Liberation clones are metrically ",
          "compatible with Arial, so line breaks will not move.")
  hit[1]
}
FIG$family <- setup_font()
message("Figure font: ", FIG$family)

# The figures are written through cairo_pdf, which resolves fonts with
# fontconfig and finds anything installed on the system. But ggplot2 also
# measures text OUTSIDE that call -- when it works out how wide an axis label is,
# how to lay out a patchwork, or where to repel a label -- and if no device is
# open it opens the one named by getOption("device"), which in a non-interactive
# session is the base pdf() device. That device resolves fonts through the
# PostScript font database, not fontconfig, so a family such as "Nimbus Sans" is
# not there and every measurement emits
#     font family 'Nimbus Sans' not found in PostScript font database
# and falls back to Helvetica's metrics.
#
# The output file is not affected -- cairo renders it correctly -- so this is not
# why a figure would come out in the wrong font. What it does affect is layout:
# widths computed from the wrong metrics move labels slightly. And fifty of these
# warnings hide the one warning that matters, which is the real cost.
#
# Registering the family in the PostScript database under Helvetica's metrics
# fixes both. Nimbus Sans and Liberation Sans are metric clones of Helvetica and
# Arial respectively, and those two are themselves metrically compatible, so the
# borrowed metrics are the right ones rather than an approximation.
register_pdf_font <- function(family) {
  if (family %in% names(grDevices::pdfFonts())) return(invisible(FALSE))
  ok <- tryCatch({
    base <- grDevices::pdfFonts()[["Helvetica"]]
    args <- setNames(list(base), family)
    do.call(grDevices::pdfFonts, args)
    do.call(grDevices::postscriptFonts, args)
    TRUE
  }, error = function(e) FALSE)
  if (ok)
    message("Registered '", family, "' in the PostScript font database ",
            "(Helvetica metrics) so text measurement stops warning.")
  else
    message("Could not register '", family, "' for text measurement; the ",
            "'not found in PostScript font database' warnings are harmless ",
            "for the output file, which is written by cairo_pdf.")
  invisible(ok)
}
register_pdf_font(FIG$family)

# --- superscripts --------------------------------------------------------------
# One rule, applied everywhere a cluster label is drawn, so that the stored label
# and the printed label cannot drift apart:
#   a trailing or embedded "+"   -> superscript plus
#   "-high" / "-low" / "-int"    -> superscript word, hyphen removed
# "Meflin+/Thy1-high" becomes "Meflin<sup>+</sup>/Thy1<sup>high</sup>".
rich <- function(x) {
  x <- as.character(x)
  x <- gsub("-(high|low|int|pos|neg)\\b", "<sup>\\1</sup>", x)
  x <- gsub("\\+", "<sup>+</sup>", x)
  x
}
use_rich <- function() isTRUE(FIG$rich_labels) && has_pkg("ggtext")
if (isTRUE(FIG$rich_labels) && !has_pkg("ggtext"))
  message("install.packages('ggtext') to draw the cluster labels with real ",
          "superscripts; until then they are drawn as stored.")

# The label function handed to a discrete scale, and the theme element that can
# render what it returns. Both fall back to plain text together -- markup drawn
# by element_text() would appear as literal <sup> tags, which is worse than no
# superscript at all.
lab_fn  <- function() if (use_rich()) rich else function(x) x

# Same labels, broken after the "/" so they occupy two short lines instead of
# one long one.
#
# CURRENTLY UNUSED. Added 2026-10-07 to stop the three subset labels colliding
# on Figure 3B's x axis and it did not work: the panel gives each subset
# 15.6 pt, and at 30 degrees the longest remaining line, "Cxcl12high", still
# needs 36 pt, so the second lines interleave much as the single lines did.
# That axis is vertical now -- see dot_plot(). This is kept because it is the
# right tool for a cramped axis that is NOT free to rotate, and because the
# split-then-mark-up order below is the part that is easy to get wrong.
#
# The split is done on the STORED name and rich() applied to each half
# afterwards, never the other way round: rich() emits "</sup>", so splitting
# the marked-up string on "/" would cut the tags in half.
#
# Labels with no "/" are returned exactly as lab_fn() would give them, which is
# why this is safe to use on an axis of gene names.
wrap_slash <- function(x) {
  f <- lab_fn()
  sep <- if (use_rich()) "<br>" else "\n"
  vapply(as.character(x), function(s) {
    if (!grepl("/", s, fixed = TRUE)) return(f(s))
    i <- regexpr("/", s, fixed = TRUE)
    paste0(f(substr(s, 1, i)), sep, f(substr(s, i + 1, nchar(s))))
  }, character(1), USE.NAMES = FALSE)
}
rich_title <- function(x) if (use_rich() && !is.null(x)) rich(x) else x
md_elem <- function(...) {
  if (use_rich()) ggtext::element_markdown(...) else element_text(...)
}

# Fixed aspect for every panel drawn on a UMAP embedding. NULL leaves ggplot to
# fill whatever space it is given, which is right for a standalone panel sized
# to fit. In a composite it is wrong twice over: the same embedding comes out
# portrait in one panel and landscape in another, so the two stop reading as the
# same kind of plot, and a stretched UMAP misrepresents the embedding geometry,
# which is the one thing a UMAP is for. fig5() and figS1() set it to 1.
#
# The cost is letterboxing: a panel whose allotted space is the wrong shape gets
# white space rather than a distorted plot. That is the right trade, but it is
# why the column split and the row heights have to be tuned together.
# Axis titles for anything drawn on a UMAP embedding. Seurat takes them from the
# reduction key, so they come out as "umap_1" or "UMAP_1" depending on how the
# object was built, and feature_plot was blanking them entirely -- hence some
# panels carrying axis names and some not. One setting, one appearance.
umap_axes <- function() {
  a <- FIG$umap_axis_titles
  if (is.null(a)) NULL else labs(x = a[1], y = a[2])
}

aspect_layer <- function() {
  a <- FIG$umap_aspect
  if (is.null(a)) NULL else theme(aspect.ratio = a)
}

theme_fig <- function(base_size = FIG$base_size) {
  theme_classic(base_size = base_size, base_family = FIG$family) +
    theme(
      axis.text        = element_text(size = base_size, colour = "black"),
      axis.title       = element_text(size = base_size + 1),
      plot.title       = element_text(size = FIG$title_size, hjust = 0),
      strip.text       = element_text(size = FIG$strip_size),
      strip.background = element_blank(),
      legend.text      = element_text(size = base_size),
      legend.title     = element_text(size = base_size),
      legend.key.size  = unit(3, "mm"))
}
theme_set(theme_fig())

# --- output -------------------------------------------------------------------
# fig() comes from section 1 of the analysis script. fig("") would return the
# path with a trailing separator, which some dir.create/file.exists combinations
# treat differently from the same path without one, so the directory is named
# without it.
# The directory this script writes to. fig() would give figures/; this must
# resolve to figures_final/ so the load message and the run banner name the
# directory the files actually land in.
FIGDIR <- out(FIG$out_subdir)

# cairo_pdf embeds the system font by name; the default pdf() device only knows
# its own base-14 set and silently maps everything else.
#
# The directory is created HERE, on every call, not once when this file is
# sourced. Creating it at load time is one state away from being wrong: delete
# the output directory in a session that already has figures.R loaded, and every
# save fails with "Cannot find directory", which is what happened after the
# unlink() that cleared the overwritten figures. Making each save responsible
# for its own directory costs nothing and cannot get out of step.
# NOTE: out(FIG$out_subdir, ...) rather than fig(). fig() is defined in the
# analysis script's section 1 and hardcodes "figures"; this script must not
# write there. Everything else about the call is unchanged.
save_plot <- function(p, name, width, height) {
  # explor_* never goes in the final directory. The exploratory fig_* functions
  # in this file call the same save_plot, and a submission folder that also
  # holds two dozen diagnostic plots is a folder someone has to sort by hand
  # before sending it. They go where the analysis script's diagnostics go.
  d <- if (startsWith(name, "explor_")) out("figures", "exploratory")
       else out(FIG$out_subdir)
  f <- file.path(d, paste0(name, ".pdf"))
  dir.create(dirname(f), showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(f, p, width = width, height = height, units = "in",
                  device = grDevices::cairo_pdf)
  # The SVG is the file to open in Illustrator (see FIG$also_svg). Same plot,
  # same size, written by a device that keeps each label as one text object and
  # names the font family instead of embedding it.
  svg_ok <- FALSE
  if (isTRUE(FIG$also_svg)) {
    if (has_pkg("svglite")) {
      fs <- file.path(d, paste0(name, ".svg"))
      # fix_text_size = FALSE is the fix for "the text box keeps its width when
      # I edit the text, and the letters stretch". svglite defaults to TRUE,
      # which writes textLength and lengthAdjust="spacingAndGlyphs" on every
      # <text> so the string occupies exactly the width R measured. Illustrator
      # honours that: change the characters and it squeezes or stretches them
      # into the old width. With FALSE the string is plain text and behaves
      # normally. The cost is that a renderer with different font metrics could
      # set the line slightly differently -- irrelevant here, because the file
      # is opened in Illustrator with the same Arial metrics it was measured
      # against.
      ggplot2::ggsave(fs, p, width = width, height = height, units = "in",
                      device = function(filename, ...)
                        svglite::svglite(filename, ..., fix_text_size = FALSE))
      svg_ok <- TRUE
    } else if (!isTRUE(.svglite_warned)) {
      message("  svglite is not installed, so only PDFs are written. ",
              "install.packages('svglite') to get Illustrator-friendly SVGs ",
              "(one text object per label, font named rather than embedded).")
      .svglite_warned <<- TRUE
    }
  }
  message("  ", name, if (svg_ok) ".pdf + .svg" else ".pdf",
          "  ", width, " x ", height, " in")
  invisible(f)
}

.svglite_warned <- FALSE

# --- raster helpers -----------------------------------------------------------
# Wraps a ggplot layer so its points become a raster while the rest stays vector.
rast <- function(layer) {
  if (FIG$raster && has_pkg("ggrastr"))
    ggrastr::rasterise(layer, dpi = FIG$raster_dpi, dev = "ragg")
  else layer
}

# Points drawn over a violin or a box. Capped per group, because drawing 1,200
# overlapping points produces a black bar, not a distribution.
# facet_col matters: the data frame is one row per cell PER GENE, so sampling by
# group alone would spread the cap across all the facets and leave a quarter of
# the points in each. Sampling by group and facet together gives every panel the
# same number.
# slice_sample() silently returns the whole group when n exceeds its size, so no
# min() is needed -- and n() cannot be used there anyway, since the n argument is
# evaluated outside the data mask.
jitter_layer <- function(df, group_col, facet_col = NULL,
                         max_n = FIG$jitter_max_n, seed = 1234) {
  d <- df
  if (!is.na(max_n)) {
    by <- c(group_col, facet_col)
    set.seed(seed)
    d <- d %>% group_by(across(all_of(by))) %>%
      slice_sample(n = max_n) %>% ungroup()
    if (nrow(d) < nrow(df))
      message("  jitter: ", nrow(d), " of ", nrow(df), " points drawn (",
              max_n, " per group", if (!is.null(facet_col)) " and panel", ")")
  }
  rast(geom_jitter(data = d, width = FIG$jitter_width, height = 0,
                   size = FIG$jitter_point, alpha = FIG$jitter_alpha,
                   colour = FIG$jitter_colour, stroke = 0))
}

# --- objects ------------------------------------------------------------------
load_obj <- function(file, what) {
  f <- obj(file)
  if (!file.exists(f))
    stop(what, " not found:\n  ", f,
         "\nRun the analysis script far enough to write it.")
  readRDS(f)
}

mouse_all <- NULL; mouse_fib <- NULL; human_all <- NULL; human_caf <- NULL

ensure <- function(what) {
  if (what == "mouse_all" && is.null(mouse_all))
    mouse_all <<- load_obj("02_allcells_annotated.rds", "annotated mouse cells")
  if (what == "mouse_fib" && is.null(mouse_fib))
    mouse_fib <<- load_obj("03c_fibroblast_subsets.rds", "mouse fibroblasts")
  if (what == "human_all" && is.null(human_all)) {
    stopifnot(!is.null(CONFIG$human_rds), file.exists(CONFIG$human_rds))
    human_all <<- readRDS(CONFIG$human_rds)
  }
  if (what == "human_caf" && is.null(human_caf)) {
    stopifnot(!is.null(CONFIG$human_caf_rds), file.exists(CONFIG$human_caf_rds))
    human_caf <<- readRDS(CONFIG$human_caf_rds)
  }
  invisible(TRUE)
}

umap_of <- function(x) {
  r <- names(x@reductions)
  hit <- r[grepl("umap", r, ignore.case = TRUE)]
  if (!length(hit)) stop("no UMAP reduction in this object")
  hit[1]
}


# =============================================================================
# UMAPs AND FEATURE PLOTS
# -----------------------------------------------------------------------------
# Seurat rasterises its own point layer when raster = TRUE, which is what keeps
# these openable in Illustrator: a UMAP of several thousand cells is several
# thousand vector circles otherwise.
# =============================================================================

# Point radius for an object of n cells, scaled from the reference and clamped.
# THE PRINTED DOT SIZE, which is the thing anyone actually looks at.
#
# scattermore draws onto a canvas of FIG$raster_px pixels and pt.size is a
# radius in canvas pixels, so a dot occupies pt/px of the panel's width
# whatever that width is. Printed diameter = 2 * pt / px * width_in.
#
# Two panels therefore match only if pt * width / px matches, NOT if pt
# matches. Figure 5A2 at 3.40 in and a Figure 5C sub-plot at 2.40 in need
# different pt values to look the same, which is why matching them by eye
# never converged.
dot_diameter_in <- function(pt, width_in, px = FIG$raster_px[1])
  2 * pt / px * width_in
pt_for_diameter <- function(diam_in, width_in, px = FIG$raster_px[1])
  diam_in * px / (2 * width_in)

pt_for <- function(n, base = FIG$umap_point, range = FIG$umap_point_range,
                   width_in = NULL, ref_in = FIG$umap_point_ref_in,
                   what = "panel") {
  p <- base * sqrt(FIG$umap_point_ref / max(n, 1))
  if (!is.null(width_in)) p <- p * width_in / ref_in
  clamped <- min(max(p, range[1]), range[2])
  # A silent clamp is how a panel ends up with dots that ignore every setting
  # you change. Say so.
  if (abs(clamped - p) > 0.05)
    message(sprintf("  NOTE: %s point radius %.1f px clamped to %.1f by ",
                    what, p, clamped),
            "the range c(", range[1], ", ", range[2], "); ",
            "widen it or change the reference width, not the point size")
  round(clamped, 1)
}

# Internal subset code ("PBF") to the label shown in the figures.
#
# The GO panel titles used to spell the labels out. They then said "Thy1-high"
# for three renamings after the object had stopped saying it, because a string
# typed in two places only agrees until one of them is edited. Reading the pair
# of columns off the object makes that impossible: whatever SUBSET_DISPLAY says
# in section 9 is what the titles say.
display_of <- function(code) {
  ensure("mouse_fib")
  md <- mouse_fib@meta.data
  if (!all(c("subset", "subset_display") %in% names(md))) return(code)
  m <- unique(data.frame(code = as.character(md$subset),
                         disp = as.character(md$subset_display)))
  out <- m$disp[match(code, m$code)]
  ifelse(is.na(out), code, out)
}

# The identities a plot should be grouped by, as a factor in plotting order.
idents_of <- function(x, group_by = NULL) {
  if (is.null(group_by)) return(Idents(x))
  v <- x[[group_by, drop = TRUE]]
  if (is.factor(v)) v else factor(v)
}

# Names or numbers? One rule, used by both the DimPlot and the feature plots, so
# the two panels of a figure cannot disagree.
label_mode <- function(ids, mode = FIG$feature_label) {
  if (!identical(mode, "auto")) return(mode)
  lv <- levels(droplevels(factor(ids)))
  if (length(lv) <= FIG$feature_label_max_n &&
      max(nchar(lv)) <= FIG$feature_label_max_ch) "name" else "number"
}

# The DimPlot is the key for any feature plot labelled with numbers. Under
# umap_key = "number" it carries the numerals on the plot and the names in the
# legend, so the key travels with the figure instead of with the caption.
umap_plot <- function(x, group_by = NULL, label = TRUE, title = NULL,
                      key = FIG$umap_key, width_in = NULL, pt = NULL) {
  # pt given explicitly overrides the size rule: used to make one panel's dots
  # print at the same size as another's.
  if (is.null(pt)) pt <- pt_for(ncol(x), width_in = width_in)
  ids <- idents_of(x, group_by)
  if (identical(key, "auto")) key <- label_mode(ids)
  message("  ", ncol(x), " cells, point radius ", pt, " px, key: ", key)

  # ADDED 2026-10-07. DimPlot's own default is scales::hue_pal(), and the
  # violin and box panels use FIG$palette, so Figure 3 was being published with
  # TWO palettes: 3A2 made Mfap5-high green and Thy1-high blue, while 3E made
  # Mfap5-high blue and Thy1-high green. A reader moving between the two panels
  # was being told the opposite thing by the colour each time. DimPlot takes
  # cols= directly, which is why this is passed as an argument rather than as a
  # second scale -- adding a scale after DimPlot replaces its own and was the
  # reason the earlier comment said not to.
  #
  # Above length(FIG$palette) levels there are not enough colours, so those
  # panels (Figure S1A, 16 clusters) keep hue_pal(). That is safe: S1A carries
  # its own 16-entry key and shares no colour meaning with any other panel.
  dim_cols <- if (nlevels(droplevels(factor(ids))) <= length(FIG$palette))
    FIG$palette[seq_len(nlevels(droplevels(factor(ids))))] else NULL

  # Seurat's own label = TRUE cannot render markup, so the on-plot labels are
  # always drawn by cluster_label_layer(): one code path, and superscripts in
  # both the name and the number modes.
  if (!identical(key, "number")) {
    p <- DimPlot(x, reduction = umap_of(x), group.by = group_by, label = FALSE,
                 pt.size = pt, cols = dim_cols, raster = FIG$raster,
                 raster.dpi = FIG$raster_px) + NoLegend()
    if (isTRUE(label))
      p <- p + cluster_label_layer(x, mode = "name", ids = ids, quiet = TRUE,
                                   size = FIG$umap_label_size)
    return(p + ggtitle(title) + theme_fig() + umap_axes() +
           theme(axis.text = element_blank(), axis.ticks = element_blank()) +
           aspect_layer())
  }

  # A throwaway metadata column whose levels are "1  Fibroblast-1", ... The
  # legend inherits the level order, so it reads 1, 2, 3, ... downwards. The
  # markup is baked into the level itself rather than applied with a scale's
  # `labels` argument, because DimPlot has already added a colour scale and a
  # second one would replace its palette.
  lv  <- levels(ids)
  key_lv <- paste0(seq_along(lv), FIG$umap_key_gap,
                   if (use_rich()) rich(lv) else lv)
  x$.umap_key <- factor(key_lv[as.integer(ids)], levels = key_lv)

  DimPlot(x, reduction = umap_of(x), group.by = ".umap_key", label = FALSE,
          pt.size = pt, cols = dim_cols, raster = FIG$raster,
          raster.dpi = FIG$raster_px) +
    cluster_label_layer(x, mode = "number", ids = ids, quiet = TRUE,
                        size = FIG$umap_label_size) +
    ggtitle(title) + theme_fig() + umap_axes() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          legend.text = md_elem(size = FIG$base_size * 0.9)) +
    guides(colour = guide_legend(override.aes = list(size = 2), ncol = 1)) +
    aspect_layer()
}

# Centroid of each cluster in the embedding, for labelling. The columns are
# named after the reduction's own axes so the labels land in the same coordinate
# space as the points.
cluster_centroids <- function(x, reduction = umap_of(x), ids = Idents(x)) {
  emb <- Embeddings(x, reduction)[, 1:2, drop = FALSE]
  ax  <- colnames(emb)                  # e.g. UMAP_1, UMAP_2
  ids <- droplevels(factor(ids))
  d <- as.data.frame(emb); d$cluster <- ids
  out <- d %>% group_by(cluster) %>%
    summarise(across(all_of(ax), median), .groups = "drop") %>%
    mutate(number = as.integer(factor(cluster, levels = levels(ids))))
  # summarise() puts the grouping variable first, so the axis columns are NOT
  # positions 1 and 2 of the result. Carrying the names through explicitly is
  # what stops the cluster column being used as the x aesthetic.
  attr(out, "axes") <- ax
  out
}

# Decides between names and numbers, then returns the layer.
# only: draw labels for these clusters alone. The numbering is still taken from
# the full set, so the numerals agree with every other panel of the same object
# -- filtering before numbering would renumber them and quietly break that.
cluster_label_layer <- function(x, mode = FIG$feature_label,
                                ids = Idents(x), quiet = FALSE,
                                size = FIG$feature_label_size, only = NULL) {
  if (identical(mode, "none")) return(NULL)
  ids  <- droplevels(factor(ids))
  mode <- label_mode(ids, mode)
  cen  <- cluster_centroids(x, ids = ids)
  # Read the attribute BEFORE touching the object: `$<-` on a tibble can drop
  # attributes it does not know about.
  ax  <- attr(cen, "axes")
  cen <- as.data.frame(cen)
  cen$label <- if (mode == "name") as.character(cen$cluster) else cen$number
  if (!is.null(only)) {
    miss <- setdiff(only, as.character(cen$cluster))
    if (length(miss))
      warning("cluster_label_layer: no such cluster: ",
              paste(miss, collapse = ", "), call. = FALSE)
    cen <- cen[as.character(cen$cluster) %in% only, , drop = FALSE]
    if (!nrow(cen)) return(NULL)
  }
  if (!quiet)
    message("  cluster labels: ", mode,
            if (mode == "number")
              paste0(" (", paste(sprintf("%d=%s", cen$number, cen$cluster),
                                 collapse = ", "), ")") else "")
  # geom_richtext renders the <sup> markup; geom_label would print the tags.
  # A name label therefore keeps its box -- markup and a halo cannot be combined
  # without a text-rendering package this script does not depend on.
  if (use_rich() && mode == "name")
    return(ggtext::geom_richtext(
      data = cen,
      aes(x = .data[[ax[1]]], y = .data[[ax[2]]], label = rich(label)),
      inherit.aes = FALSE, size = size, family = FIG$family,
      label.colour = NA, label.padding = unit(rep(0.4, 4), "mm"),
      fill = alpha("white", 0.65), colour = "black"))

  if (!identical(FIG$label_style, "halo"))
    return(geom_label(data = cen,
                      aes(x = .data[[ax[1]]], y = .data[[ax[2]]], label = label),
                      inherit.aes = FALSE, size = size,
                      family = FIG$family, linewidth = 0,
                      label.padding = unit(0.4, "mm"),
                      fill = alpha("white", 0.65), colour = "black"))

  # Halo, drawn by shadowtext when it is installed. That package renders the
  # outline in device units, so the halo is the same thickness on a 3.5-inch
  # UMAP and on a 0.9-inch feature plot in the same figure -- which is the whole
  # point, since the glyph is a fixed number of points on both.
  if (has_pkg("shadowtext"))
    return(shadowtext::geom_shadowtext(
      data = cen,
      aes(x = .data[[ax[1]]], y = .data[[ax[2]]], label = label),
      inherit.aes = FALSE, size = size, family = FIG$family,
      colour = "black", bg.colour = "white", bg.r = FIG$label_halo_r))

  # Fallback: eight white copies offset in DATA units, then the black glyph.
  # Approximate, because a data-unit offset scales with the panel while the text
  # does not, so the halo is relatively thinner on a small panel. Still better
  # than a box, which covers exactly the cells the label is pointing at.
  if (!isTRUE(.halo_warned)) {
    message("  shadowtext is not installed; drawing an approximate halo. ",
            "install.packages('shadowtext') for a device-unit outline.")
    .halo_warned <<- TRUE
  }
  rng <- apply(cen[, ax], 2, function(v) diff(range(v)))
  d   <- FIG$label_halo_frac * mean(rng)
  ang <- seq(0, 2 * pi, length.out = 9)[-9]
  halo <- lapply(ang, function(a) {
    cc <- cen
    cc[[ax[1]]] <- cc[[ax[1]]] + d * cos(a)
    cc[[ax[2]]] <- cc[[ax[2]]] + d * sin(a)
    geom_text(data = cc,
              aes(x = .data[[ax[1]]], y = .data[[ax[2]]], label = label),
              inherit.aes = FALSE, size = size, family = FIG$family,
              colour = "white")
  })
  c(halo,
    list(geom_text(data = cen,
                   aes(x = .data[[ax[1]]], y = .data[[ax[2]]], label = label),
                   inherit.aes = FALSE, size = size, family = FIG$family,
                   colour = "black")))
}
.halo_warned <- FALSE

# Per-gene colour limits, and the numbers behind them. Returns a data frame with
# one row per gene so the choice can be checked rather than trusted: pct is the
# percentage of cells with any detection, and hi is where the colour saturates.
# A gene whose pct is a fraction of a per cent is a gene whose panel will look
# empty however the scale is set -- that is a fact about the data, and it should
# be read off this table rather than inferred from the picture.
expr_cutoffs <- function(x, features, q = FIG$feature_hi_quantile,
                         min_cells = FIG$feature_min_express) {
  m <- FetchData(x, vars = features, layer = "data")
  out <- lapply(features, function(f) {
    v  <- m[[f]]
    nz <- v[v > 0]
    hi <- if (!length(nz)) NA_real_
          else if (length(nz) < min_cells) max(nz)
          else as.numeric(stats::quantile(nz, q))
    data.frame(gene = f, n_pos = length(nz),
               pct = 100 * length(nz) / length(v),
               hi = hi, max = if (length(nz)) max(nz) else 0)
  })
  do.call(rbind, out)
}

# group_by must match the companion DimPlot's, or the numerals on the two panels
# mean different clusters. Passing it explicitly at every call site is what stops
# that: Idents() and the celltype column are not guaranteed to agree after a
# re-run.
feature_plot <- function(x, features, ncol = 3, label = FIG$feature_label,
                         width_in = NULL,
                         group_by = NULL) {
  f <- intersect(features, rownames(x))
  missing <- setdiff(features, f)
  if (length(missing)) message("  absent: ", paste(missing, collapse = ", "))
  # width_in here is the width of the WHOLE panel; each sub-plot gets a
  # ncol-th of it, and that is what sets the dot size.
  sub_w <- if (is.null(width_in)) NULL else width_in / ncol
  pt <- pt_for(ncol(x), FIG$feature_point, FIG$feature_point_range,
               width_in = sub_w, ref_in = FIG$feature_point_ref_in,
               what = "feature plot")
  message("  ", ncol(x), " cells, point radius ", pt, " px",
          if (!is.null(sub_w)) sprintf("  (sub-plot %.2f in wide)", sub_w))

  # order = TRUE draws the expressing cells last, on top of the grey ones. It is
  # not optional here: with a rare population the positive cells are painted over
  # by the negative ones otherwise.
  if (identical(FIG$feature_scale, "expressing")) {
    cut <- expr_cutoffs(x, f)
    message("  colour scale from expressing cells only:")
    for (i in seq_len(nrow(cut)))
      message(sprintf("    %-8s %6d pos (%5.2f%%)  saturates at %.2f of max %.2f",
                      cut$gene[i], cut$n_pos[i], cut$pct[i],
                      cut$hi[i], cut$max[i]))
    lo_cut <- rep(0, nrow(cut))
    hi_cut <- ifelse(is.na(cut$hi), 1, cut$hi)
  } else {
    lo_cut <- FIG$feature_min_cutoff
    hi_cut <- FIG$feature_max_cutoff
  }

  p <- FeaturePlot(x, features = f, reduction = umap_of(x), ncol = ncol,
                   cols = FIG$feature_cols, order = TRUE,
                   min.cutoff = lo_cut, max.cutoff = hi_cut,
                   pt.size = pt, raster = FIG$raster,
                   raster.dpi = FIG$raster_px)
  # cluster_label_layer now returns a LIST of layers in halo mode (eight white
  # copies plus the black one). patchwork's `&` takes one element at a time, so
  # add them one by one rather than handing it the list.
  lab <- cluster_label_layer(x, label, ids = idents_of(x, group_by))
  if (!is.null(lab)) {
    ll <- if (inherits(lab, "gg") || inherits(lab, "Layer")) list(lab) else lab
    for (l in ll) p <- p & l
  }
  p <- p &
    theme_fig() &
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          plot.title = element_text(size = FIG$title_size, face = "italic"))
  ax <- umap_axes()
  p <- if (is.null(ax)) p & theme(axis.title = element_blank()) else p & ax
  al <- aspect_layer()
  if (!is.null(al)) p <- p & al
  p
}

fig_umaps <- function() {
  message("UMAPs and feature plots")
  ensure("mouse_all"); ensure("mouse_fib")

  # Wider than the fibroblast panels: the numbered key adds a legend column.
  save_plot(umap_plot(mouse_all, group_by = "celltype"),
            "explor_umap_allcells", 6.0, 3.8)
  save_plot(feature_plot(mouse_all, c("Islr", "Thy1", "Col1a1",
                                      "Pdgfra", "Pdgfrb", "Acta2"),
                         group_by = "celltype"),
            "explor_featureplots_allcells", 6.6, 4.4)

  grp <- if ("subset_display" %in% colnames(mouse_fib@meta.data))
    "subset_display" else "subset"
  # key = "name" is pinned rather than left to "auto": the labels are now long
  # enough to sit near the auto threshold, and crossing it would silently swap
  # this panel from on-plot names to numbers plus a legend.
  save_plot(umap_plot(mouse_fib, group_by = grp, key = "name"),
            "explor_umap_fibroblasts", 3.6, 3.2)
  save_plot(feature_plot(mouse_fib, c("Islr", "Thy1", "Pi16",
                                      "Cxcl12", "Col1a1", "Acta2"),
                         group_by = grp),
            "explor_featureplots_fibroblasts", 6.6, 4.4)

  if (!is.null(CONFIG$human_rds) && file.exists(CONFIG$human_rds)) {
    ensure("human_all")
    save_plot(umap_plot(human_all), "explor_umap_human_allcells", 6.0, 3.8)
    save_plot(feature_plot(human_all, c("ISLR", "THY1", "COL1A1",
                                        "PDGFRA", "PDGFRB", "ACTA2")),
              "explor_featureplots_human_allcells", 6.6, 4.4)
  }
  if (!is.null(CONFIG$human_caf_rds) && file.exists(CONFIG$human_caf_rds)) {
    ensure("human_caf")
    save_plot(umap_plot(human_caf), "explor_human_CAF_umap", 3.6, 3.2)
    save_plot(feature_plot(human_caf, c("ISLR", "THY1", "PI16",
                                        "CXCL12", "COL1A1", "ACTA2")),
              "explor_featureplots_human_CAF", 6.6, 4.4)
  }
  invisible(TRUE)
}


# =============================================================================
# DOT PLOTS
# =============================================================================

# scale = TRUE is Seurat's default: each gene is z-scored ACROSS the identities,
# so the colour shows where a gene is highest relative to itself. That is the
# right choice for a marker panel, where the genes have wildly different
# absolute levels and the question is which cell type each one marks.
#
# It is the wrong choice when the genes are being compared with one another, and
# it becomes actively misleading when there are few identities: with three
# groups a z-score has three possible values, so a gene detected in 0.3 % of
# every group is painted with the same reds as a gene detected in 48 % of one
# group and 0.6 % of the others. Pass scale = FALSE there; the colour is then
# mean log-normalised expression, comparable between genes, and an undetected
# gene is a pale dot instead of a false positive.
# flip = TRUE puts the genes on the vertical axis and the identities on the
# horizontal one, which turns a very wide short panel into a tall narrow one.
# That is the shape that lets the dot plot sit beside the enrichment panels
# instead of spanning the figure on its own row.
dot_plot <- function(x, features, order_levels = NULL, dot_scale = 3,
                     scale = TRUE, flip = FALSE,
                     legend = c("right", "top")) {
  legend <- match.arg(legend)
  f <- intersect(features, rownames(x))
  missing <- setdiff(features, f)
  if (length(missing)) message("  absent: ", paste(missing, collapse = ", "))
  if (!is.null(order_levels)) {
    keep <- intersect(order_levels, levels(Idents(x)))
    # An order that matches nothing used to be applied anyway, turning every
    # identity into NA and producing a plot that was wrong rather than absent.
    # This is how FigS2B was drawn against numbered clusters after the cell-type
    # column it expected was missing from the object.
    if (!length(keep))
      stop("none of the requested identity levels exist in this object. ",
           "Asked for: ", paste(utils::head(order_levels, 5), collapse = ", "),
           if (length(order_levels) > 5) ", ..." else "",
           " | present: ", paste(utils::head(levels(Idents(x)), 8),
                                 collapse = ", "))
    if (length(keep) < length(order_levels))
      message("  not in this object, dropped from the order: ",
              paste(setdiff(order_levels, keep), collapse = ", "))
    Idents(x) <- factor(as.character(Idents(x)), levels = keep)
  }
  lv <- levels(Idents(x))
  # DotPlot warns when scaling across only a few identities. With four CAF
  # clusters the warning is expected and harmless.
  p <- suppressWarnings(DotPlot(x, features = f, dot.scale = dot_scale,
                                scale = scale,
                                cols = c("grey90", FIG$palette[1]))) +
    theme_fig() + theme(axis.title = element_blank())

  # legend = "top" puts both guides in one horizontal band above the panel
  # instead of a column to the right of it. ADDED 2026-10-08 for Figure S1C,
  # where the two guides were taking about 90 pt -- a sixth -- of a 540 pt
  # page, squeezing 44 gene columns into roughly 390 pt. Overhead, the band
  # costs height rather than width and the dots get the whole page.
  #
  # The band costs about 0.35 in of height, so a panel that moves its guides
  # up should be given that height back (Figure S1C went 3.5 -> 3.9 in) or the
  # rows close up on each other.
  #
  # justification "left" so the band starts at the page edge rather than
  # floating in the middle above the panel.
  if (identical(legend, "top"))
    p <- p +
      guides(
        colour = guide_colourbar(title.position = "left", title.vjust = 1,
                                 barwidth  = unit(1.1, "in"),
                                 barheight = unit(0.09, "in"), order = 1),
        size   = guide_legend(title.position = "left", nrow = 1, order = 2)) +
      theme(legend.position      = "top",
            legend.box           = "horizontal",
            legend.direction     = "horizontal",
            legend.justification = "left",
            legend.box.spacing   = unit(3, "pt"),
            legend.margin        = margin(t = 0, b = 0))

  if (!flip)
    return(p +
      scale_y_discrete(limits = rev(lv), labels = lab_fn()) +
      theme(axis.text.y = md_elem(size = FIG$base_size, colour = "black"),
            axis.text.x = element_text(angle = FIG$dot_angle, hjust = 1,
                                       vjust = 1, face = "italic")))
  # Flipped. DotPlot draws features along x in the order given, so after
  # coord_flip() they would read bottom-to-top; rev() puts the first gene at the
  # top, where a reader looks for it. The identities are NOT reversed here --
  # coord_flip already does that for the y axis.
  # CHANGED 2026-10-01: the identities were rev(lv) here. With coord_flip that
  # put them on the horizontal axis in the reverse of the factor order, so this
  # panel read CxF, PPF, PBF while every other panel read PBF, PPF, CxF. The
  # order now comes from the factor levels, which are set once in the analysis
  # script, so all of Fig. 3 reads the same way.
  p +
    scale_x_discrete(limits = rev(f)) +
    scale_y_discrete(limits = lv, labels = lab_fn()) +
    coord_flip() +
    # THE SUBSET LABELS ARE VERTICAL, and the reason is arithmetic rather than
    # taste. Measured off the 2026-10-07 PDF, this panel gives each of the
    # three subsets 15.6 pt of axis, and "Meflin+/Cxcl12high" sets about 70 pt
    # long. A label at angle t occupies 70*cos(t) pt of the axis, so it only
    # fits inside 15.6 pt once cos(t) < 0.22, that is beyond 77 degrees. At 30
    # degrees it needs 61 pt and the three labels were drawn on top of one
    # another; breaking them after the "/" was tried on 2026-10-07 and did not
    # help either, because the longest remaining line, "Cxcl12high", still
    # needs 36 pt and the second lines then interleave as well.
    #
    # At 90 degrees the label occupies its LINE HEIGHT on the axis, about
    # 10 pt, which fits inside 15.6 pt with room to spare. It costs height --
    # the longest label is 70 pt of the 360 pt page -- and that is the cheapest
    # of the three currencies here, the others being the gene column (the
    # panel's whole point) and the dots.
    #
    # Two conditions to re-check if this panel is ever resized: a wider panel
    # may let the labels back down to 30 degrees, and a SECOND line at
    # 90 degrees would need about 19 pt and would start colliding again.
    #
    # The left plot.margin added earlier is gone with the angle that needed it:
    # vertical labels have no leftward overhang, and those 10 pt go back to the
    # dots, which only have about 47 pt to work with.
    theme(axis.text.y = element_text(size = FIG$base_size, colour = "black",
                                     face = "italic"),
          axis.text.x = md_elem(size = FIG$base_size, colour = "black",
                                angle = 90, hjust = 1, vjust = 0.5))
}

MOUSE_DOT_GENES <- c(
  "Islr", "Col1a1", "Col1a2", "Pdgfra", "Pdgfrb", "Dcn", "Lum",
  "Thy1", "Pi16", "Mfap5", "Scara5",
  "Lrat", "Reln", "Lama1", "Rgs5", "Notch3", "Cspg4",
  "Acta2", "Tagln", "Cnn1", "Myh11", "Des", "Pln",
  "Pecam1", "Cdh5", "Stab2", "Clec4g", "Vwf", "Tek",
  "Sema3g", "Ly6c1", "Lyve1", "Prox1",
  "Msln", "Upk3b", "Lrrn4", "Krt19", "Krt7",
  # ADDED 2026-09-24: Epcam, Sox9 and Cftr. Without them this panel shows that
  # the keratin-positive cluster is Msln+ but not that it is Epcam-, and the
  # second half is what distinguishes mesothelium from biliary epithelium --
  # mesothelium expresses Krt19 and Krt7 too. The panel is cited in the Results
  # for the annotation of that cluster, so it has to carry the negative as well
  # as the positive.
  "Epcam", "Sox9", "Cftr",
  "Alb", "Fabp1", "Apoa5", "Rgn", "Ptprc", "Cd3e")

HUMAN_DOT_GENES <- c(
  "ISLR", "COL1A1", "COL1A2", "PDGFRA", "PDGFRB", "DCN", "LUM",
  "THY1", "PI16", "MFAP5", "SCARA5",
  "RGS5", "NOTCH3", "CSPG4",
  "ACTA2", "TAGLN", "CNN1", "MYH11", "DES",
  "PECAM1", "CDH5", "VWF", "TEK", "LYVE1", "PROX1", "PDPN",
  "EPCAM", "KRT19", "KRT7", "CDH1",
  "PLP1", "S100B", "SOX10", "NGFR",
  "PTPRC", "CD3E", "NKG7", "CD68", "LYZ", "C1QA",
  "MS4A1", "CD79A", "TPSAB1", "CPA3")

# Markers of the three fibroblast subsets, for the panel that Results paragraph
# 2 cites. Grouped by subset and then by what the genes are, so the blocks read
# in the same order as the sentence.
#
# The last two blocks are the basement-membrane collagens and the secreted
# morphogen antagonists. They belong HERE, in a marker panel, and not in a
# module score or a separate hypothesis panel: a marker dot plot makes no claim
# beyond "these genes are restricted to this subset", which is what the sentence
# says and what the differential-expression table supports. The moment the same
# genes are collapsed into a score and given a P value, the circularity starts
# -- the genes were chosen because they are markers of that subset.
# Every gene here has to earn its column against one test: does a claim in the
# Results depend on the reader seeing it? Twenty-five columns was a list of
# everything interesting rather than a panel.
#
# DROPPED, with reasons:
#   Dkk2, Hhip   the morphogen-antagonist reading is not developed anywhere in
#                the paper now: the panel that scored it was circular and is
#                gone, the module score panel (Fig. 5E) no longer holds it,
#                and the Discussion
#                does not connect it to Meflin. Four columns for a claim nobody
#                makes is four columns spent.
#   Clec3b       97.1 % against 85.2 %; a gene present in nearly every cell of
#                every subset shows the reader nothing.
#   Gja4         below the expression floor in the subset identification, and
#                Notch3 already carries the perivascular block.
#   Fbn2, Ntrk2  markers of the Thy1-high subset with nothing resting on them.
#
# KEPT, against the same test:
#   Notum, Wif1  these are the two most subset-restricted genes in the whole
#                dataset -- 43.7 % and 48.6 % against 0.6-1.6 % elsewhere, 36-
#                and 23-fold, far above Thy1's 2.6-fold. Dropping the sharpest
#                markers of a subset from its marker panel because the pathway
#                story was abandoned would be the wrong correction: they are
#                shown as markers, and the text lists them as markers, with no
#                claim about Wnt attached.
#   Col4a5/6,    the peribiliary position rests on these; collagen IV and VII are
#     Col7a1     basement-membrane specific, which is an argument independent of
#                the clustering.
#   Pi16         no longer a name, but now a finding: shared by the Thy1-high and
#                Mfap5+ subsets and absent from Cxcl12-high, which is the pattern
#                the human data repeat.
#   Pdgfrb       added for continuity with the in situ panels in Figures 1-4,
#                which use it.
# =============================================================================
# THE PUBLISHED FIBROBLAST FRAMEWORK (Buechler et al., Nature 2021, 593:575-579,
# doi:10.1038/s41586-021-03549-5)
# -----------------------------------------------------------------------------
# That study integrated ~230,000 fibroblasts from 17 tissues and named two
# UNIVERSAL subtypes, Pi16+ and Col15a1+, from which tissue-specialised and
# activated states derive. Every subtype is named after its defining gene.
#
# PROVENANCE, AND WHAT THIS PANEL IS NOT. The paper's full per-subtype marker
# lists are in its supplementary tables, which are not open access; they could
# not be retrieved here. So this is NOT a module score of the published
# signatures -- it is the expression of the genes the published subtypes are
# NAMED AFTER. That is weaker as a similarity measure and stronger as a
# citation: every gene here is verifiable from the paper's figures, and nothing
# is attributed to a list I could not read.
#
# To upgrade it to a proper module score, download the supplementary marker
# tables (institutional access) and hand them over; the scoring machinery in
# section 11 takes gene lists directly.
#
# CHECK THE ROSTER BEFORE IT GOES IN A LEGEND. I am confident about Pi16,
# Col15a1, Cxcl12, Ccl19, Lrrc15, Comp, Npnt, Hhip, Adamdec1 and Bmp4; less so
# that Fbln1 and Coch belong and that none is missing. Confirm against Fig. 1
# and Fig. 3 of the paper.
BUECHLER_GENES <- c(
  # universal
  "Pi16", "Col15a1",
  # tissue-specialised, steady state
  "Npnt", "Hhip", "Comp", "Bmp4", "Cxcl12", "Adamdec1", "Coch", "Fbln1",
  # activated / perturbed
  "Lrrc15", "Ccl19")

FIB_DOT_GENES <- c(
  "Islr", "Col1a1", "Pdgfra", "Pdgfrb",
  # Dkk2 and Hhip ADDED 2026-10-05, beside Notum and Wif1, which were already
  # here. The Results now name the Thy1-high subset's programme as epithelial
  # morphogenesis carried by morphogen antagonists; with only two of the four
  # drawn, a reader cannot see the block. Detection in the subset against the
  # other two: Notum 43.7/1.1, Wif1 48.6/1.4, Dkk2 42.4/1.0, Hhip 36.0/2.4 --
  # Dkk2 ranks 33rd and Hhip 50th of this subset's 527 depth-robust markers by
  # that difference, alongside Notum (30th) and Wif1 (21st), so they belong to
  # the same tier and are not a weaker pair added to make a point.
  #
  # NOT a module score and NOT a new gene panel. Figure S1F was deleted on
  # 2026-09-16 because its Wnt/Hedgehog panel had been assembled from this
  # subset's own top markers, so showing the genes were restricted to it tested
  # nothing. These four are individual marker genes in the subset marker dot
  # plot, which is what that panel already is; no score is computed and no
  # enrichment is claimed from them.
  "Thy1", "Ptn", "Notum", "Wif1", "Dkk2", "Hhip",
  "Col4a5", "Col4a6", "Col7a1",
  "Pi16", "Mfap5", "Scara5", "Pcolce2", "Fbln2",
  # Gja4 added 2026-09-24: 33.0 % vs 2.6 % detection, 27.3-fold on linear means,
  # the largest separation of any gene tested for this subset. It is wrong for a
  # violin (median 0) and right for a dot plot, where the detection rate is the
  # dot size and a gene in a third of the cells still reads.
  # Colec11 ADDED 2026-10-05. It is this subset's LARGEST separator of all 779
  # depth-robust markers (93.6 % against 28.8 %, log2FC 4.06) and was absent
  # while two of the other subsets' top separators were shown; under the
  # selection rule now stated in the Methods its omission was the anomaly. It
  # is also an innate-immune lectin, which is consistent with -- but not the
  # reason for -- its inclusion: the rule is the detection-rate difference.
  "Cxcl12", "Colec11", "Notch3", "Heyl", "Fbln5", "Gja4"
  # ADDED then REMOVED 2026-09-29/30. Mki67 and Top2a were added here to
  # document an ABSENCE (1.9 / 0.6 / 3.9 % and 0.6 / 0.3 / 1.5 % detection) and
  # taken out again the next day. Two reasons, and the second is the fatal one:
  #   1. This panel is the subset MARKER dot plot. Neither gene is a marker of
  #      any subset, so their presence changes what the panel is.
  #   2. Colour is z-scored within each gene, so a gene detected in two per cent
  #      of cells still gets the full colour range and the Cxcl12-high subset
  #      reads as "high" for Mki67 -- the opposite of the claim. A panel that
  #      needs its legend to say "ignore the colour for these two rows" is
  #      mis-designed, not merely in need of a caption.
  # The absence is already carried where it does work: the Figure 3D legend
  # states Mki67 <= 3.9 % and Top2a <= 1.5 %, next to the GO term
  # "mesenchymal cell proliferation" that it exists to disarm.
  #
  # Acta2 was considered for the same block and rejected for the same scaling
  # reason: 70.7 / 42.2 / 48.8 % detection would draw a large dot and, after
  # scaling, make the Thy1-high subset read as Acta2-high. The absence of a
  # myofibroblastic population is a comparison against smooth muscle, not within
  # the fibroblasts, and FigS1C already draws it: Acta2, Tagln, Cnn1, Myh11, Des
  # and Pln across the fibroblast clusters AND VSM-1, VSM-2, Mural and Pericyte.
  # Cite that panel instead.
  )

# Genes for the violin panels, mouse and human.
#
# These two lists are deliberately the same list in two species, because the two
# panels are read against each other: the cross-species argument is that THY1
# and the adventitial markers, which are co-expressed in one murine subset,
# separate between two human clusters. A reader cannot check that if the mouse
# panel shows Thy1 and Cxcl12 while the human panel shows THY1 and MFAP5.
#
# What each block is doing, and what it is NOT doing:
#   Islr           the paper's claim that all subsets are Meflin-positive
#   Thy1/Pi16/     the genes the subsets are NAMED after. Showing that the
#     Cxcl12         Thy1-high subset is Thy1-high is true by construction, so
#                    these panels are there to show the SHAPE of the
#                    distribution -- whether the marker is bimodal or graded,
#                    and whether the subsets overlap. They are not evidence of
#                    anything and must not carry a P value in the text.
#   Mfap5/Scara5   the adventitial ECM programme. This is the informative part
#                    of the panel: it is what splits between the two human
#                    clusters and it was not used to define anything.
#   Col1a1/Acta2   fibroblast identity and the absence of a myofibroblastic
#                    state, which is a separate claim in the Results.
#   Col4a5/Col4a6/ the basement-membrane collagens. Added because the
#     Col7a1       peribiliary claim rests on them and the panel otherwise
#                    contains nothing that supports it: Thy1 is the name of the
#                    subset, so it cannot also be the evidence for what the
#                    subset is. Collagen IV and VII are basement-membrane
#                    specific, which is an anatomical argument independent of
#                    the clustering. They are markers of the subset, so they
#                    carry the same caveat as any marker -- fine in a display,
#                    not as a scored programme with a P value.
# Twelve genes, three per subset, in a 4 x 3 grid.
#
# The original eleven were unevenly spread: four for the Thy1-high subset (311
# cells) and one for the Cxcl12-high subset, which is 1,165 cells, 65 % of the
# fibroblasts. Col4a6 was dropped -- COL4A5 and COL4A6 share a bidirectional
# promoter on Xq22 and are co-expressed, so the second panel duplicated the
# first. Col7a1 is a different collagen and stays.
#
# THE TWO ADDED Cxcl12-high GENES, AND WHY NOT CHEMOKINES.
# The obvious choice was to follow the GO result and show chemokines, since
# chemotaxis is that subset's strongest term (58 genes, FE 3.24, adjusted
# P 1.8e-12). The detection rates say no. Measured in mouse_fib by subset:
#
#     gene     CxF det   next    linear sep   CxF median
#     Pf4        40.3 %   9.3 %        6.7x        0
#     Cx3cl1     55.5    27.3          2.2         0.64
#     Ccl19      64.8    34.7          2.3         0.83
#     Cxcl5      36.7    18.2          2.1         0
#     Cxcl3      25.4    12.9          3.0         0
#     Cxcl1      97.6    97.1          1.4        40.66
#     Ccl7       92.7    87.3          --         15.45
#
# Only Pf4 clears 30 % detection with a threefold separation, and even its
# median is zero. A violin is a picture of a distribution; genes detected in a
# minority of cells draw as a flat line with a few dots, and three of those
# beside a strong GO panel would make the enrichment look like nothing. The
# chemotaxis term is real because 58 genes move together, not because any one
# chemokine marks the subset -- which is itself worth a sentence in the text.
# Note also that Cxcl1 and Ccl7 are expressed comparably in ALL THREE subsets,
# so neither belongs in a list of subset-specific chemokines.
#
# What is shown instead:
#   Ednrb   63.3 % vs 7.5 %, 11.3x, median 0.96 vs 0 -- the single best
#           discriminating gene for this subset, better than Notch3 (9.9x).
#   Isg15   93.6 % vs 65.0 %, 3.0x, median 7.68 vs 1.38 -- separates by MEDIAN,
#           which is what a violin can show, and it is part of the type I
#           interferon module that the GO panel reports, so the two panels tell
#           one story without either being a chemokine near its detection floor.
#
# Isg15 rather than Iigp1, and NOT both: they are both interferon-stimulated
# genes, so two panels would say "this subset is interferon-high" twice -- the
# same duplication that Col4a5 and Col4a6 were removed for. Isg15 is the more
# visible of the two (median 7.68 against 2.39). Iigp1 is the swap if a less
# ubiquitous ISG is wanted: FIG-free, just edit the list below.
#
# Ednrb and Iigp1 are also less loaded than Notch3, which a reader would query
# as a mural marker. The perivascular reading is not claimed in the Results, so
# the figure should not push it.
VLN_GENES_MOUSE <- c("Islr", "Col1a1",                        # pan-fibroblast
                     "Thy1", "Col4a5", "Col7a1",              # Thy1-high
                     "Pi16", "Mfap5", "Scara5",               # Mfap5+
                     "Cxcl12", "Ednrb", "Isg15",              # Cxcl12-high
                     "Acta2")                                 # myofibroblast

VLN_GENES_HUMAN <- toupper(VLN_GENES_MOUSE)

# ADDED 2026-10-07 -- Figure S1D (S1E until the 2026-10-08 renumber), gene by gene.
#
# WHAT IT IS FOR. Figure 3B carries 24 genes in a flipped dot plot: each row is
# a few millimetres high, and its colour is z-scored WITHIN each gene, so the
# panel shows which subset is highest but not by how much. Every gene the text
# leans on is in there, and none of them is legible as a quantity. This panel
# redraws the ones a statement depends on, on a real axis, with the cells
# shown. The supplement is where a reader goes to check a specific gene, so the
# overlap with 3B is the point rather than a cost.
#
# THE SELECTION RULE: every gene on which a sentence in the single-cell Results
# or in the Discussion's fibroblast paragraph depends, MINUS the six that
# already have a feature plot of their own in Figure 3C (Islr, Col1a1, Pi16,
# Thy1, Mfap5, Cxcl12) and the pan-fibroblast references no claim rests on
# (Pdgfra, Pdgfrb). Fifteen genes remain, and none of them is currently
# legible anywhere except as a row of the 24-gene dot plot.
#
# Pi16 is excluded under that rule even though it heads one of the signatures;
# the legend points at Figure 3C so the reader is not left looking for it.
#
# It also closes a gap: Dpt, Ly6c1, Sema3c and Ackr3 appeared in NO figure at
# all, while the Results say "Its individual components ... were also most
# highly expressed in that subset". Figure 3E shows the module SCORE, an
# aggregate, which cannot support a per-gene sentence -- and the same paragraph
# says Pi16 alone does not separate the Mfap5-high subset from the Thy1-high
# one, so a reader is entitled to ask whether the score rests on one or two
# genes.
#
# WHY VIOLINS AND NOT A SECOND DOT PLOT. A second dot plot would repeat 3B's
# problem: colour z-scored per gene draws Dpt (99.0 % against 97.8 %) across
# the same range as Sema3c (25.0 % against 4.2 %). That is the argument that
# removed Mki67, Top2a and Acta2 from Figure 3B. Free-y violins with the cells
# overlaid show the distributions as they are. For the near-binary genes
# (Notum, Wif1, Dkk2, Hhip) the violin is mostly a mass at zero and the
# information is in the points and the median bar -- the legend says so.
# CHANGED 2026-10-08. The blocks were named after the SUBSET the genes mark
# ("Thy1-high", "Mfap5-high") for two of them and after the published
# SIGNATURE they belong to ("Pi16+ universal") for the other two -- two
# different kinds of reason on one axis. The subset names were also the
# uninformative half: the colour key already says which subset is which, so a
# strip reading "Thy1-high" told the reader nothing the violin had not.
#
# The blocks now say WHAT THE GENE IS. For the Thy1-high markers that is a
# function, and it lines the panel up one-for-one with the Discussion sentence
# "distinguished by basement-membrane collagens and by secreted antagonists of
# Wnt and Hedgehog signalling (Notum, Wif1, Dkk2, Hhip)". The two signature
# blocks keep their names, because for those genes membership of a published
# signature IS what they are here for. "Secreted" is left to the legend: it
# applies to all four antagonists and would not fit on a strip.
#
# SCARA5 REMOVED, 2026-10-08, 15 genes -> 14. It was the one gene whose block
# could not be named this way. Its three possible meanings each failed:
#   - function (ferritin receptor, Li et al. Dev Cell 2009) -- this manuscript
#     makes no claim about Scara5's function, and a functional strip label
#     would insinuate one;
#   - "adventitial", the pre-specified signature it belongs to -- that
#     programme was reclassified data_derived on 2026-09-29 and dropped from
#     Figure 3E because its genes are markers of the subset it scores highest
#     in, so naming it here would quietly restore a retired claim;
#   - "Mfap5-high marker" -- the uninformative kind of label this change is
#     removing.
# Nothing is orphaned by the removal: the only sentence that names Scara5 is
# the Discussion's "tracing driven by a subset-restricted gene such as Scara5
# or Mfap5", and its restriction is already visible in Figure 3B, where dot
# SIZE is percent expressed (44.8 % against 1.1 % in the other two subsets).
# Mfap5, the other gene that sentence names, is a feature plot in Figure 3C.
# INTERFERON BLOCK ADDED 2026-10-08, 14 genes -> 18, and the grid 5 -> 6
# columns. The Results say "the strongest enrichment in the Meflin+/Cxcl12high
# subset was an interferon-response programme" and that sentence had NO
# gene-level figure: only the GO term in Fig. 3D and the NES in Fig. 3E. Same
# gap, and same remedy, as Dpt/Ly6c1/Sema3c/Ackr3.
#
# HOW THE FOUR WERE CHOSEN, and why not the obvious ones. Picking the ISGs with
# the biggest fold change out of the Cxcl12-high marker list is the mistake
# that got Figure S1F deleted on 2026-09-16: genes taken from a subset's marker
# list are high in that subset by construction. So the starting set is the
# intersection of the two interferon results the manuscript actually reports --
# the genes of GO:0035456 "response to interferon-beta", the term DRAWN in
# Figure 3D, and the leading edge of HALLMARK_INTERFERON_ALPHA_RESPONSE, the
# set drawn in Figure 3E. Six genes: Ifit3, Gbp3, Bst2, Ifitm3, Irf1, Ifitm1.
#
# Two are then dropped on DISPLAY grounds, which the legend has to state:
#   Ifitm3  100 % of Cxcl12-high cells against 99.8 % of the others -- a violin
#           cannot show a difference that is not in the detection rate, and a
#           flat facet reads as evidence against the claim it illustrates
#   Ifitm1  also a depth-robust marker of the Thy1-high subset (97.4 % vs
#           84.8 %), so it is not a subset-restricted gene at all
#
# BE HONEST ABOUT THIS: two clauses applied to six genes landing on exactly the
# four that fill a 6 x 3 grid is a coincidence worth declaring rather than
# hiding. Each clause is defensible on its own and neither uses the effect size
# in this dataset, but the rule was written knowing the answer.
#
# What the row will look like, from Table S2 (detection in the subset vs the
# other two): Ifit3 42.4/13.6, Gbp3 34.9/18.7, Bst2 93.1/78.0, Irf1 82.8/64.3.
# Only Ifit3 separates strongly. That is the honest picture of a programme
# called by GSEA across many genes rather than by a few switch-like ones, and
# it is why the block is in the supplement and not in Figure 3.
NAMED_GENE_SETS <- list(
  `Wnt antagonist`      = c("Notum", "Wif1", "Dkk2"),
  `Hedgehog antagonist` = c("Hhip"),
  `Basement membrane`   = c("Col4a5", "Col7a1"),
  `Interferon response` = c("Ifit3", "Gbp3", "Bst2", "Irf1"),
  `Pi16+ universal`     = c("Dpt", "Ly6c1", "Sema3c", "Ackr3"),
  `Col15a1+ universal`  = c("Col15a1", "Hspg2", "Vwa1", "Cd34"))

# The order above is the order the Results introduce the findings, and with six
# columns it also puts the whole Thy1-high block -- the three antagonists and
# the two basement-membrane collagens, six genes -- on the first row by itself.
# The other blocks have four genes each and cannot align with a six-column
# grid; that is what the strip labels are for.

# How each block name is set. The names above are the plain-text form and are
# what the strip shows when ggtext is absent; the values are the marked-up form
# for element_markdown(). Only the two signature blocks differ: Pi16 and
# Col15a1 are gene symbols and belong in italic with a superscript plus, the
# same as everywhere else in the figure. A block with no entry here is drawn as
# its name.
NAMED_GENE_BLOCK_LABELS <- c(
  `Pi16+ universal`    = "*Pi16*<sup>+</sup> universal",
  `Col15a1+ universal` = "*Col15a1*<sup>+</sup> universal")

# --- how compact Figure S1D is drawn -----------------------------------------
# ADDED 2026-10-07. The first drawn version was 6.8 x 5.0 in and did not fit on
# the Figure S1 page. Almost all of that height went to one thing: the three
# subset names, angled at 30 degrees, repeated under all fifteen facets. They
# carry no information that one legend cannot carry, so NAMED_GENE_COMPACT moves them
# into a single horizontal key under the panel. The facets then no longer have
# to be as wide as the longest angled label either, so the panel loses width as
# well as height.
#
#   NAMED_GENE_COMPACT = TRUE,  NAMED_GENE_BLOCK_STRIP = TRUE   5.8 x 3.55 in   (the panel)
#   NAMED_GENE_COMPACT = TRUE,  NAMED_GENE_BLOCK_STRIP = FALSE  4.6 x 3.15 in
#   NAMED_GENE_COMPACT = FALSE                          6.8 x 5.00 in   (as first drawn)
#
# NAMED_GENE_BLOCK_STRIP keeps the second strip line naming what each gene is. It is
# the reason the panel exists -- the contrast between the blocks -- so it is on
# by default. Turning it off buys another 1.2 in of width and 0.4 in of height,
# and the grouping then HAS to be stated in the figure legend instead, because
# the facet order alone does not show where one block ends: with five columns
# the blocks do not line up with the rows.
#
# The 2026-10-08 removal of Scara5 leaves 14 genes, so the third row of the
# 5-column grid holds four facets and one empty cell. That is deliberate: it
# keeps the Col15a1+ block whole on its own row, which five columns of 15 did
# not do.
NAMED_GENE_COMPACT    <- TRUE
NAMED_GENE_BLOCK_STRIP <- TRUE
# Columns in the grid. Six for the eighteen genes: three full rows, with the
# whole Thy1-high block on the first. IF THE GENE LIST CHANGES, change this and
# PANEL_SIZE$named_gene_vln together -- the width is about 1.1 in per column
# plus 0.3 for the y axis, and named_gene_size() cannot work that out for you.
NAMED_GENE_NCOL       <- 6

# Candidate genes for naming the subsets, checked with name_check() below.
NAME_CANDIDATES_MOUSE <- c("Islr", "Thy1", "Pi16", "Cxcl12", "Mfap5", "Scara5",
                           "Clec3b", "Dpt", "Pcolce2", "Notch3", "Notum")

# CAF state markers, for checking whether the human clusters deserve the names
# they carry. The lists are the published definitions, taken from the papers
# that introduced the terms rather than assembled from this dataset:
#   iCAF / myCAF   Ohlund et al., J Exp Med 2017; Elyada et al., Cancer Discov
#                  2019 (the CXCL/IL and the ACTA2/TAGLN programmes respectively)
#   apCAF          Elyada et al. 2019 (MHC class II)
# A cluster called iCAF that is negative for the whole iCAF block is misnamed,
# and the name then propagates into every sentence that cites the cluster.
CAF_STATE_GENES <- list(
  iCAF      = c("IL6", "IL11", "LIF", "CXCL1", "CXCL2", "CXCL8", "CXCL12",
                "HAS1", "HAS2", "C3", "CFD", "DPT", "PDGFRA"),
  myCAF     = c("ACTA2", "TAGLN", "POSTN", "TNC", "THBS2", "CTHRC1", "TGFB1"),
  apCAF     = c("HLA-DRA", "HLA-DRB1", "CD74", "SLPI"),
  universal = c("PI16", "ISLR", "THY1", "MFAP5", "SCARA5", "CLEC3B"),
  # Not a CAF state. Added because the unbiased markers of the Meflin_ssl
  # cluster -- CDH19 in 79 % of its cells against 1 % elsewhere, DCX, COL28A1,
  # CLDN1, P2RY14 -- are Schwann-cell and perineurial genes, not fibroblast
  # genes. A cluster of 48 cells carrying that profile has to be excluded as a
  # neural-lineage population before it can be reported as a fibroblast subset,
  # and it is the cluster the cross-species comparison depends on.
  schwann   = c("PLP1", "S100B", "SOX10", "MPZ", "NGFR", "CDH19", "MBP",
                "CLDN1", "COL28A1", "DCX", "PMP22", "L1CAM"),
  # Lineage sentinels, the same idea used on the mouse object. The human CAF
  # object was subset from the whole tissue by a label; nothing since has
  # checked that everything in it is a fibroblast. These blocks ask that
  # directly, and they are what the naming argument should have started from.
  fibroblast  = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FBLN1",
                  "PDGFRA", "VIM"),
  mural       = c("RGS5", "PDGFRB", "NOTCH3", "KCNJ8", "HIGD1B", "ACTA2",
                  "MYH11", "DES", "CSPG4"),
  perineurial = c("CLDN1", "SLC2A1", "ITGB4", "EMP2", "GJA1", "MUC1",
                  "KRT18", "ANGPTL7"),
  other       = c("EPCAM", "KRT19", "PECAM1", "PTPRC", "MSLN", "UPK3B"))

# Unbiased markers of each human CAF cluster, and how well each cluster is
# separated from its nearest neighbour.
#
# The blocks in caf_check() answer "does this cluster fit a name we already
# have". This answers the prior question: what does this cluster actually
# express that the others do not. It is the analysis the original annotation
# should have rested on, and running it is cheaper than arguing about a label.
#
# n_sep counts genes passing the thresholds against the NEAREST other cluster,
# not against all others pooled. A cluster that separates from the pool but not
# from its neighbour is a subdivision of that neighbour, and no name will fix
# that -- it is a clustering-resolution question.
caf_markers <- function(n = 15, min_pct = 0.25, logfc = 0.5) {
  ensure("human_caf")
  x <- human_caf
  lv <- levels(droplevels(Idents(x)))
  message("clusters: ", paste(sprintf("%s (n=%d)", lv, table(Idents(x))[lv]),
                              collapse = ", "))

  m <- Seurat::FindAllMarkers(x, only.pos = TRUE, min.pct = min_pct,
                              logfc.threshold = logfc, verbose = FALSE)
  for (cl in lv) {
    all_cl <- m[m$cluster == cl, ]
    # Ranking by fold change alone puts genes with an adjusted P of 1.0 at the
    # top of the list, which is how a cluster comes to look well defined by
    # markers that are not significant. Significant genes first, then fold
    # change within them.
    sig <- all_cl[all_cl$p_val_adj < 0.05, ]
    top <- sig[order(-sig$avg_log2FC), ][seq_len(min(n, nrow(sig))), ]
    message("\n=== ", cl, ": ", nrow(all_cl), " markers vs all others, ",
            nrow(sig), " with adjusted P < 0.05; top ", nrow(top),
            " of those by log2FC")
    for (i in seq_len(nrow(top)))
      message(sprintf("    %-12s log2FC %5.2f  pct %4.2f / %4.2f  p.adj %8.1e",
                      top$gene[i], top$avg_log2FC[i], top$pct.1[i],
                      top$pct.2[i], top$p_val_adj[i]))
  }

  message("\npairwise separation (genes with |log2FC| > ", logfc,
          " and adjusted P < 0.05):")
  pairs <- utils::combn(lv, 2, simplify = FALSE)
  for (p in pairs) {
    d <- tryCatch(Seurat::FindMarkers(x, ident.1 = p[1], ident.2 = p[2],
                                      min.pct = min_pct,
                                      logfc.threshold = logfc,
                                      verbose = FALSE),
                  error = function(e) NULL)
    k <- if (is.null(d)) NA_integer_ else sum(d$p_val_adj < 0.05)
    message(sprintf("  %-12s vs %-12s %5s genes", p[1], p[2],
                    ifelse(is.na(k), "error", k)))
  }
  invisible(m)
}

# Cross-tabulation of two markers, per cluster and over the whole object.
#
# This is the readout an immunostain gives -- double positive, single positive
# either way, double negative -- computed from the sequencing so the two can be
# compared directly instead of argued about. A subset-level summary cannot
# answer it: "ISLR is detected in 77 % and THY1 in 89 %" is compatible with
# every cell being double positive and with the two sets barely overlapping.
#
# Dropout matters here and is not a detail. ISLR sits at 0.8-2.3 mean normalised
# counts in this object against 21-49 for COL1A1, and at that level a large
# share of expressing cells return zero. The double-negative and single-negative
# fractions below are therefore UPPER bounds, and an immunostain -- protein,
# whole cell, amplified -- will legitimately show more positive cells than this
# table does. It is the ORDERING and the spatial pattern that transfer between
# the two methods, not the percentages.
coexpr_check <- function(x = "human_caf", g1 = "ISLR", g2 = "THY1",
                         group_by = NULL) {
  x <- as_obj(x)
  if (is.null(group_by)) {
    group_by <- intersect(c("subset_display", "subset"), colnames(x@meta.data))
    group_by <- if (length(group_by)) group_by[1] else NULL
  }
  ids <- droplevels(idents_of(x, group_by))
  f <- c(g1, g2)
  if (!all(f %in% rownames(x))) {
    message("absent: ", paste(setdiff(f, rownames(x)), collapse = ", "))
    return(invisible(NULL))
  }
  d <- FetchData(x, vars = f, layer = "data")
  a <- d[[g1]] > 0
  b <- d[[g2]] > 0

  message(sprintf("co-expression of %s and %s, %% of cells in each cluster",
                  g1, g2))
  message(sprintf("  %-22s %6s %8s %8s %8s %7s",
                  "cluster", "n", paste0(g1, "+", g2, "+"),
                  paste0(g1, "+", g2, "-"), paste0(g1, "-", g2, "+"), "both-"))
  for (cl in c(levels(ids), "ALL")) {
    k <- if (identical(cl, "ALL")) rep(TRUE, length(a)) else ids == cl
    n <- sum(k)
    message(sprintf("  %-22s %6d %7.1f%% %7.1f%% %7.1f%% %6.1f%%",
                    cl, n,
                    100 * sum(a[k] & b[k]) / n, 100 * sum(a[k] & !b[k]) / n,
                    100 * sum(!a[k] & b[k]) / n, 100 * sum(!a[k] & !b[k]) / n))
  }
  invisible(data.frame(cluster = ids, g1 = a, g2 = b))
}

# Which patients each cluster came from.
#
# 364 CAFs from 7 cases is about 52 per tumour, and at that recovery a cluster
# of 48 cells can easily be one patient's cells. A population present in one
# case is a case report; a population present in six is a finding. Nothing else
# in this audit distinguishes the two, so this is the check that decides how
# much weight each cluster can carry -- and it applies to the perineurial
# cluster and to the myCAF/iCAF split equally.
#
# Prints the candidate sample columns first, because the column holding the
# patient identifier is named differently in every object.
sample_check <- function(x = "human_caf", sample_col = NULL, group_by = NULL) {
  x <- as_obj(x)
  md <- x@meta.data
  cand <- names(md)[vapply(md, function(v)
    (is.character(v) || is.factor(v)) &&
      dplyr::between(dplyr::n_distinct(v), 2, 30), logical(1))]
  message("categorical columns with 2-30 levels: ",
          if (length(cand)) paste(cand, collapse = ", ") else "(none)")
  if (is.null(sample_col)) {
    hit <- grep("orig.ident|sample|patient|case|donor|subject|library",
                cand, ignore.case = TRUE, value = TRUE)
    sample_col <- if (length(hit)) hit[1] else NULL
  }
  if (is.null(group_by)) {
    group_by <- intersect(c("subset_display", "subset"), names(md))
    group_by <- if (length(group_by)) group_by[1] else NULL
  }
  ids <- droplevels(idents_of(x, group_by))

  # No sample column: fall back to the barcode suffix. Seurat's merge() appends
  # _1 ... _n (or -1 ... -n from CellRanger's aggregation) to keep barcodes
  # unique, so the suffix IS the library of origin even when nothing recorded it
  # in the metadata. This is a last resort, not a substitute for the real
  # annotation: the suffix says which library a cell came from, and whether one
  # library is one patient is something only the person who built the object
  # knows. The mapping is printed so it can be checked against the sample sheet.
  if (is.null(sample_col)) {
    suf <- sub("^.*[-_]", "", colnames(x))
    if (length(unique(suf)) > 1 && length(unique(suf)) <= 30) {
      message("no sample column; using the barcode suffix instead — ",
              length(unique(suf)), " libraries: ",
              paste(sort(unique(suf)), collapse = ", "))
      message("CHECK that one suffix is one patient before using this.")
      s <- factor(suf, levels = sort(unique(suf)))
    } else {
      message("no sample column and the barcode suffix is not usable; ",
              "pass sample_col = \"<name>\"")
      return(invisible(NULL))
    }
  } else {
    s <- factor(md[[sample_col]])
  }
  tb  <- table(ids, s)

  message("\ncells per cluster and ",
          if (is.null(sample_col)) "library (barcode suffix)" else sample_col,
          ":")
  message("  ", sprintf("%-22s", ""),
          paste(sprintf("%8s", colnames(tb)), collapse = ""), "   total")
  for (i in seq_len(nrow(tb)))
    message("  ", sprintf("%-22s", rownames(tb)[i]),
            paste(sprintf("%8d", tb[i, ]), collapse = ""),
            sprintf("%8d", sum(tb[i, ])))

  message("\nhow concentrated each cluster is:")
  for (i in seq_len(nrow(tb))) {
    v <- tb[i, ]
    message(sprintf("  %-22s %2d of %2d samples contribute; largest is %.0f%%",
                    rownames(tb)[i], sum(v > 0), length(v),
                    100 * max(v) / sum(v)))
  }
  invisible(tb)
}

# The stored clusterings, cross-tabulated against the labels in use.
#
# The human object carries RNA_snn_res.0.1 through 1.5 from when it was built.
# That makes the myCAF/iCAF question answerable without re-clustering anything:
# if the two labels fall into a single cluster at a lower resolution, then the
# object's own clustering already says they are one population, and merging
# them is restoring a split that a higher resolution introduced rather than
# overruling the analysis.
resolution_check <- function(x = "human_caf", group_by = NULL) {
  x <- as_obj(x)
  md  <- x@meta.data
  res <- grep("^RNA_snn_res\\.", names(md), value = TRUE)
  res <- res[order(as.numeric(sub("^RNA_snn_res\\.", "", res)))]
  if (!length(res)) {
    message("no RNA_snn_res.* columns in this object")
    return(invisible(NULL))
  }
  if (is.null(group_by)) {
    group_by <- intersect(c("subset_display", "subset"), names(md))
    group_by <- if (length(group_by)) group_by[1] else NULL
  }
  ids <- droplevels(idents_of(x, group_by))
  for (r in res) {
    tb <- table(ids, factor(md[[r]]))
    tb <- tb[, colSums(tb) > 0, drop = FALSE]
    message("\n", r, " (", ncol(tb), " clusters)")
    message("  ", sprintf("%-22s", ""),
            paste(sprintf("%7s", colnames(tb)), collapse = ""))
    for (i in seq_len(nrow(tb)))
      message("  ", sprintf("%-22s", rownames(tb)[i]),
              paste(sprintf("%7d", tb[i, ]), collapse = ""))
    # Two labels that always land in the same column at this resolution are one
    # population as far as this clustering is concerned.
    home <- apply(tb, 1, function(v) colnames(tb)[which.max(v)])
    dup  <- home[duplicated(home) | duplicated(home, fromLast = TRUE)]
    if (length(dup))
      message("  labels sharing a cluster here: ",
              paste(sprintf("%s->%s", names(dup), dup), collapse = ", "))
  }
  invisible(TRUE)
}

# What the whole-tissue object called the cells of a CAF cluster.
#
# The `annotations` column inside the CAF object turned out to be a copy of the
# CAF labels, so it cannot say what these cells were before the fibroblast
# subset was taken. The whole-tissue object can, by barcode.
#
# Barcodes are matched exactly first. If the two objects were built with
# different suffixes the exact match returns nothing, so the prefix (the 16-base
# barcode) is tried as well -- that is ambiguous across libraries and is
# reported as such rather than used silently.
trace_cells <- function(cluster, from = NULL, to = NULL, to_col = "annotations") {
  if (is.null(from)) { ensure("human_caf"); from <- human_caf }
  if (is.null(to))   { ensure("human_all"); to   <- human_all }
  if (!to_col %in% colnames(to@meta.data)) {
    message(to_col, " not in the whole-tissue object; columns: ",
            paste(names(to@meta.data), collapse = ", "))
    return(invisible(NULL))
  }
  ids <- colnames(from)[droplevels(Idents(from)) == cluster]
  if (!length(ids)) {
    message("no cells with identity ", cluster)
    return(invisible(NULL))
  }
  hit <- intersect(ids, colnames(to))
  message(cluster, ": ", length(ids), " cells, ", length(hit),
          " matched in the whole-tissue object by exact barcode")
  if (!length(hit)) {
    pf <- sub("[-_].*$", "", ids)
    pt <- sub("[-_].*$", "", colnames(to))
    k  <- pt %in% pf
    message("  exact match failed; matching on the barcode prefix instead ",
            "gives ", sum(k), " cells for ", length(unique(pf)),
            " distinct prefixes -- AMBIGUOUS across libraries, treat as a hint")
    if (!sum(k)) return(invisible(NULL))
    tb <- table(to@meta.data[[to_col]][k])
  } else {
    tb <- table(to@meta.data[[to_col]][match(hit, colnames(to))])
  }
  tb <- sort(tb[tb > 0], decreasing = TRUE)
  for (i in seq_along(tb))
    message(sprintf("    %-28s %4d  (%.0f%%)", names(tb)[i], tb[i],
                    100 * tb[i] / sum(tb)))
  invisible(tb)
}

# Runs the naming check on the human CAF object, one published block at a time.
# Blocks are kept separate because the question is per block: is the cluster
# called iCAF the one that scores highest on the iCAF genes, and is it doing so
# across the block or on one gene?
caf_check <- function(blocks = names(CAF_STATE_GENES)) {
  ensure("human_caf")
  out <- list()
  for (b in blocks) {
    message("\n=== ", b, " markers ===")
    out[[b]] <- name_check(human_caf, CAF_STATE_GENES[[b]])
  }
  invisible(out)
}

# The human CAF panel, mirroring FIB_DOT_GENES block for block so that the two
# species can be read against each other. Deriving it from the mouse list rather
# than writing it out again is deliberate: two hand-maintained lists drift, and
# the whole point of the panel is that they are the same genes.
#
# Two blocks are then APPENDED, because the human object contains something the
# mouse object does not. The mouse tissue is normal and has no activated
# fibroblasts; the human tissue is tumour and two of its four clusters are
# activated states. Without these genes, myCAFs and iCAFs would appear in the
# panel with nothing that defines them, which invites the reader to conclude
# they are negative for everything rather than that they were not asked about.
FIB_DOT_GENES_HUMAN <- c(
  toupper(FIB_DOT_GENES),
  "ACTA2", "TAGLN", "POSTN",          # myofibroblastic
  "IL6", "CXCL1", "CXCL8", "HAS1")    # inflammatory

# Cluster order for that panel: the gene blocks run from the Thy1-high
# phenotype through the adventitial one to the activated states, so the clusters
# are ordered to match and the panel reads along its diagonal. Names that do not
# exist in the object are dropped by dot_plot(), so this is safe to leave as is
# if the labels are renamed.
HUMAN_CAF_ORDER <- c("Thy1_ssl", "Meflin_ssl", "myCAF", "iCAF")

HUMAN_CELLTYPE_ORDER <- c(
  "Fibroblast", "Smooth_muscle_cell", "Endothelial_cell",
  "Lymphatic_endothelial_cell", "Epithelial_cell", "Schwann_cell",
  "T_cell", "NK_cell", "T_reg", "Macrophage", "B_cell", "Mast_cell")

fig_dotplots <- function() {
  message("Dot plots")
  ensure("mouse_all")
  Idents(mouse_all) <- mouse_all$celltype
  save_plot(dot_plot(mouse_all, MOUSE_DOT_GENES), "explor_dotplot_allcells",
            8.0, 3.2)

  # The marker panel for the three subsets. This is the display Results
  # paragraph 2 cites for the basement-membrane collagens and the morphogen
  # antagonists, alongside the differential-expression table.
  ensure("mouse_fib")
  x <- mouse_fib
  Idents(x) <- idents_of(x, if ("subset_display" %in% colnames(x@meta.data))
    "subset_display" else "subset")
  # Flipped, to match the composite: genes down the side, identities across.
  save_plot(dot_plot(x, FIB_DOT_GENES, dot_scale = 3, flip = TRUE),
            "explor_dotplot_fibroblast_subsets", 2.6, 4.2)

  if (!is.null(CONFIG$human_rds) && file.exists(CONFIG$human_rds)) {
    ensure("human_all")
    grp <- CONFIG$human_group_by
    if (!is.null(grp) && grp %in% colnames(human_all@meta.data)) {
      Idents(human_all) <- human_all@meta.data[[grp]]
    } else if (!is.null(grp)) {
      # The fallback to Idents() is usually right -- this object carries the
      # cell-type labels there rather than in a metadata column -- but it was
      # silent, so there was no way to tell a correct fallback from a wrong one.
      # Printing the identities it fell back to makes that a one-line check.
      message("  CONFIG$human_group_by = '", grp, "' is not a column here; ",
              "using Idents(): ",
              paste(head(levels(droplevels(Idents(human_all))), 12),
                    collapse = ", "))
    }
    save_plot(dot_plot(human_all, HUMAN_DOT_GENES, HUMAN_CELLTYPE_ORDER),
              "explor_dotplot_human_allcells", 8.0, 3.0)
  }
  if (!is.null(CONFIG$human_caf_rds) && file.exists(CONFIG$human_caf_rds)) {
    ensure("human_caf")
    save_plot(dot_plot(human_caf, FIB_DOT_GENES_HUMAN, HUMAN_CAF_ORDER,
                       dot_scale = 4),
              "explor_human_CAF_dotplot", 6.8, 2.2)
  }
  invisible(TRUE)
}


# =============================================================================
# VIOLINS AND BOX PLOTS
# -----------------------------------------------------------------------------
# Built with ggplot rather than VlnPlot so that the point layer can be capped,
# made transparent and rasterised independently of the violin outline.
# =============================================================================

# x_axis and legend were added 2026-10-07 for the named-gene panel. Both default to the
# behaviour every existing caller already had, so nothing else in this file
# changes. x_axis = FALSE drops the per-facet group labels and legend = TRUE
# puts the same three names in one horizontal key underneath: in a many-facet
# panel the labels are repeated once per facet and the key is drawn once, which
# is where the height goes. Use them together -- x_axis = FALSE on its own
# leaves the violins unlabelled.
violin_panel <- function(df, group_col, value_col, facet_col = NULL,
                         ylab = NULL, points = TRUE, free_y = TRUE,
                         ncol = NULL, x_axis = TRUE, legend = FALSE) {
  p <- ggplot(df, aes(.data[[group_col]], .data[[value_col]],
                      fill = .data[[group_col]])) +
    geom_violin(scale = "width", linewidth = FIG$violin_linewidth,
                colour = "grey20", alpha = 0.85, trim = TRUE)
  if (points) p <- p + jitter_layer(df, group_col, facet_col)
  p <- p +
    stat_summary(fun = median, geom = "crossbar", width = 0.5,
                 linewidth = 0.25, colour = "black") +
    # labels = lab_fn() applies the same rich-text rendering to the key that
    # the x axis gets; with legend = FALSE (every other caller) it is inert.
    scale_fill_manual(values = FIG$palette, labels = lab_fn()) +
    scale_x_discrete(labels = lab_fn()) +
    labs(x = NULL, y = ylab) + theme_fig() +
    theme(legend.position = "none",
          axis.text.x = md_elem(angle = 30, hjust = 1,
                                size = FIG$base_size, colour = "black"),
          # see dot_plot(): the angled labels overhang the left edge of the
          # page. Harmless when x_axis = FALSE; the compact branch of
          # fig_named_gene_violins() sets its own margin after this anyway.
          plot.margin = margin(5.5, 5.5, 5.5, 16))
  if (!x_axis)
    p <- p + theme(axis.text.x = element_blank(),
                   axis.ticks.x = element_blank())
  if (legend)
    # override.aes undoes the violins' alpha so the key swatches are solid;
    # without it they read as a paler version of the fills they stand for.
    p <- p + labs(fill = NULL) +
      guides(fill = guide_legend(nrow = 1, override.aes = list(alpha = 1))) +
      theme(legend.position = "bottom",
            legend.direction = "horizontal",
            legend.margin = margin(t = 0, b = 0),
            legend.box.spacing = unit(4, "pt"),
            # theme_fig() sets legend.text with element_text(), which would
            # print the superscript markup lab_fn() produces as literal tags.
            legend.text = md_elem(size = FIG$base_size, colour = "black"))
  if (!is.null(facet_col))
    p <- p + facet_wrap(vars(.data[[facet_col]]), ncol = ncol,
                        scales = if (free_y) "free_y" else "fixed")
  p
}

# ncol and aspect were added 2026-09-29 for Fig. 5E. aspect = 1 forces every
# facet's plotting region square regardless of how the outer canvas is sized,
# which is the only way to guarantee square facets when the y scales are free
# and the strip heights differ; the cost is that the saved page may be
# letterboxed, which is harmless because the panels are assembled by hand.
# Both default to NULL, so every other caller is unchanged.
box_panel <- function(df, group_col, value_col, facet_col = NULL, ylab = NULL,
                      ncol = NULL, aspect = NULL) {
  p <- ggplot(df, aes(.data[[group_col]], .data[[value_col]],
                      fill = .data[[group_col]]))
  # points first, box over them at partial opacity, so nothing is hidden
  p <- p + jitter_layer(df, group_col, facet_col) +
    geom_boxplot(outlier.shape = NA, alpha = FIG$box_alpha,
                 linewidth = FIG$box_linewidth, colour = "grey20") +
    scale_fill_manual(values = FIG$palette) +
    scale_x_discrete(labels = lab_fn()) +
    labs(x = NULL, y = ylab) + theme_fig() +
    theme(legend.position = "none",
          axis.text.x = md_elem(angle = 30, hjust = 1,
                                size = FIG$base_size, colour = "black"),
          # see dot_plot(): the angled labels overhang the left edge. In the
          # 2026-10-07 Fig. 3E the "Me" of the first label was 2.7 pt outside
          # the page box.
          plot.margin = margin(5.5, 5.5, 5.5, 16))
  if (!is.null(facet_col))
    p <- p + facet_wrap(vars(.data[[facet_col]]), scales = "free_y", ncol = ncol)
  if (!is.null(aspect)) p <- p + theme(aspect.ratio = aspect)
  p
}

# group_col is a metadata column name, or NULL to use the object's identities.
expr_long <- function(x, features, group_col = NULL) {
  f <- intersect(features, rownames(x))
  missing <- setdiff(features, f)
  # A gene silently dropped from a panel is worse than an error: the figure
  # still looks complete. Report what is missing and what the object does have
  # under a similar name, which is usually a case or a suffix difference.
  if (length(missing)) {
    message("  ABSENT from this object: ", paste(missing, collapse = ", "))
    for (g in missing) {
      hit <- grep(paste0("^", substr(g, 1, 4)), rownames(x),
                  ignore.case = TRUE, value = TRUE)
      if (length(hit))
        message("    ", g, " -> similar names present: ",
                paste(head(hit, 8), collapse = ", "))
    }
  }
  g <- if (is.null(group_col) || !group_col %in% colnames(x@meta.data))
    as.character(Idents(x)) else as.character(x@meta.data[[group_col]])
  m <- tryCatch(GetAssayData(x, layer = "data"),
                error = function(e) GetAssayData(x, slot = "data"))
  m <- as.matrix(m[f, , drop = FALSE])
  as_tibble(t(m)) %>%
    mutate(group = g) %>%
    pivot_longer(-group, names_to = "gene", values_to = "expression") %>%
    mutate(gene = factor(gene, levels = f))
}

# =============================================================================
# NAMING EVIDENCE
# -----------------------------------------------------------------------------
# Not a figure. Prints what is needed to decide what to call a subset, which is
# not the same question as which gene is significant.
#
# A name should be the gene that SEPARATES the subset from its neighbours, and
# it should be expressed in most of its cells, or the name misdescribes the
# majority of the population it labels. Those two properties are what this
# prints: pct is the percentage of cells in the subset with any detection, and
# the separation columns are the ratio of the top subset to the next one, in
# detection and in mean expression.
#
# A gene can be the subset's most significant marker and still be a poor name:
# significance rises with cell number, and a gene detected in 30 % of a subset
# against 5 % elsewhere is highly significant while leaving 70 % of the cells
# unaccounted for by the name.
# =============================================================================

# x may be a Seurat object OR the name of one ("mouse_all", "mouse_caf", ...),
# which is loaded if it is not in the session yet.
#
# Naming it as a string is the safer call and the reason this exists. The
# objects start as NULL and are filled in by ensure(), so name_check(mouse_all)
# BEFORE anything has loaded mouse_all passes NULL, which the old default then
# silently replaced with the fibroblast object -- a table of the wrong cells
# under the right heading. Passing "mouse_all" cannot do that.
as_obj <- function(x) {
  if (is.character(x) && length(x) == 1) { ensure(x); return(get(x, envir = globalenv())) }
  if (is.null(x))
    stop("the object is NULL. It has not been loaded yet: pass its name as a ",
         "string, e.g. name_check(\"mouse_all\"), or call ensure() first.")
  x
}

name_check <- function(x = "mouse_fib", genes = NAME_CANDIDATES_MOUSE,
                       group_by = NULL) {
  x <- as_obj(x)
  # The mouse object carries its labels in a metadata column; the human CAF
  # object, which came from elsewhere, carries them in Idents() and has no
  # "subset" column at all. Falling through to Idents() rather than naming a
  # column that may not exist is what lets this be called on either.
  if (is.null(group_by)) {
    group_by <- intersect(c("subset_display", "subset"), colnames(x@meta.data))
    group_by <- if (length(group_by)) group_by[1] else NULL
  }
  ids <- droplevels(idents_of(x, group_by))
  message("grouping by ", if (is.null(group_by)) "Idents()" else group_by,
          ": ", paste(levels(ids), collapse = ", "))
  f   <- intersect(genes, rownames(x))
  if (length(setdiff(genes, f)))
    message("absent: ", paste(setdiff(genes, f), collapse = ", "))
  d <- FetchData(x, vars = f, layer = "data")

  out <- lapply(f, function(g) {
    v   <- d[[g]]
    pct <- tapply(v > 0, ids, mean) * 100
    mu  <- tapply(v, ids, mean)              # mean of log-normalised values
    lin <- tapply(expm1(v), ids, mean)       # mean of normalised counts
    o   <- order(-pct)
    ol  <- order(-lin)
    data.frame(gene = g, top = names(pct)[o[1]],
               pct_top = pct[o[1]], pct_next = pct[o[2]],
               sep_pct = pct[o[1]] / max(pct[o[2]], 0.1),
               mean_top = mu[o[1]],
               sep_mean = mu[o[1]] / max(mu[order(-mu)][2], 1e-3),
               # The linear-scale mean is what avg_expr() reports elsewhere in
               # this project. It is dominated by the few highest cells, so it
               # can rank the subsets differently from the log-scale mean and
               # from detection. Where the three disagree, the gene does not
               # separate the subsets and should not be used to name one.
               lin_top = names(lin)[ol[1]],
               lin_val = lin[ol[1]],
               sep_lin = lin[ol[1]] / max(lin[ol][2], 1e-3))
  })
  out <- do.call(rbind, out)
  out$agree <- out$top == out$lin_top
  out <- out[order(out$top, -out$sep_pct), ]

  message("naming evidence:")
  message(sprintf("  %-9s %-12s %7s %8s %7s %8s %8s  %-12s %8s",
                  "gene", "highest pct", "pct", "pct next", "sep", "mean",
                  "sep mean", "highest lin", "sep lin"))
  for (i in seq_len(nrow(out)))
    message(sprintf("  %-9s %-12s %6.1f%% %7.1f%% %6.1fx %8.2f %7.1fx  %-12s %7.1fx%s",
                    out$gene[i], out$top[i], out$pct_top[i], out$pct_next[i],
                    out$sep_pct[i], out$mean_top[i], out$sep_mean[i],
                    out$lin_top[i], out$sep_lin[i],
                    if (!out$agree[i]) "  <- METRICS DISAGREE" else ""))
  # "The metrics disagree" is not one finding, it is two, and the earlier version
  # of this function conflated them. Detection and level disagree for a gene that
  # separates the subsets by neither -- and also for a gene detected in nearly
  # every cell of every subset, where detection is saturated and carries no
  # information at all while the level still differs several-fold. The first is a
  # reason not to use the gene; the second is a reason to use its level and
  # ignore its detection rate. Islr is the second kind and was being reported as
  # the first.
  min_det <- apply(vapply(f, function(g) tapply(d[[g]] > 0, ids, mean) * 100,
                          numeric(length(levels(ids)))), 2, min)
  out$min_det <- min_det[out$gene]
  out$verdict <- ifelse(
    out$min_det > 80,
    "detected in every subset: judge on level, not detection",
    ifelse(!out$agree & pmax(out$sep_pct, out$sep_lin) < 1.5,
           "separates nothing: do not use",
           ifelse(!out$agree, "metrics disagree: check the medians", "")))
  flagged <- out[nzchar(out$verdict), ]
  if (nrow(flagged)) {
    message("")
    for (i in seq_len(nrow(flagged)))
      message("  ", sprintf("%-9s", flagged$gene[i]), flagged$verdict[i])
  }

  # The summary above collapses the comparison to top-versus-next, which hides
  # WHICH subset the runner-up is. For a naming decision that is the thing you
  # need: a gene detected in 50 %, 36 % and 5 % of the three subsets is a very
  # different proposition from 50 %, 48 % and 5 %. Print every cell.
  lv <- levels(ids)
  full_pct <- vapply(f, function(g) tapply(d[[g]] > 0, ids, mean) * 100,
                     numeric(length(lv)))
  full_lin <- vapply(f, function(g) tapply(expm1(d[[g]]), ids, mean),
                     numeric(length(lv)))
  show_matrix <- function(m, title, fmt) {
    message("\n  ", title)
    message("    ", sprintf("%-10s", ""),
            paste(sprintf("%12s", colnames(m)), collapse = ""))
    for (i in seq_len(nrow(m)))
      message("    ", sprintf("%-10s", rownames(m)[i]),
              paste(sprintf(fmt, m[i, ]), collapse = ""))
  }
  # The median says whether a higher mean is the whole population moving up or a
  # minority of cells pulling the average. Only the first justifies calling a
  # subset "X-high": if the medians are equal and the means are not, the name
  # describes a tail, not the cells it labels.
  full_med <- vapply(f, function(g) tapply(expm1(d[[g]]), ids, median),
                     numeric(length(lv)))
  show_matrix(full_pct, "detection rate (% of cells), all subsets", "%11.1f%%")
  show_matrix(full_lin, "mean normalised counts, all subsets", "%12.2f")
  show_matrix(full_med, "MEDIAN normalised counts, all subsets", "%12.2f")
  invisible(list(summary = out, pct = full_pct, lin = full_lin,
                 median = full_med))
}

fig_violins <- function() {
  message("Violin and box plots")
  ensure("mouse_fib")
  grp <- if ("subset_display" %in% colnames(mouse_fib@meta.data))
    "subset_display" else "subset"
  lv <- levels(mouse_fib@meta.data[[grp]])

  d <- expr_long(mouse_fib, VLN_GENES_MOUSE, grp) %>%
    mutate(group = factor(group, levels = lv))
  save_plot(violin_panel(d, "group", "expression", "gene",
                         ylab = "Normalised expression") +
              theme(strip.text = element_text(face = "italic")),
            "explor_violin_subset_markers", 6.4, 5.6)

  sc <- mouse_fib@meta.data %>%
    select(matches("^programme_[0-9]+$")) %>% as_tibble()
  if (ncol(sc)) {
    # AddModuleScore names the columns programme_1 ... programme_n in the order
    # of the list it was given, and that order is written to
    # Fig5_07_programme_definitions.csv together with the data_derived flag.
    # Reading it back does two things: the facets get real names instead of
    # "programme_3", and the programmes built from this dataset's own markers are
    # kept out of the main figure. A score whose genes are the cluster's markers
    # cannot be evidence about that cluster, so it goes to its own supplementary
    # panel, labelled as descriptive.
    defs <- tab("Fig5_07_programme_definitions.csv")
    nm <- NULL
    if (file.exists(defs)) {
      d <- read.csv(defs, stringsAsFactors = FALSE)
      nm <- setNames(d$programme, paste0("programme_", seq_len(nrow(d))))
      derived <- setNames(as.logical(d$data_derived),
                          paste0("programme_", seq_len(nrow(d))))
    } else {
      message("  Fig5_07_programme_definitions.csv not found; ",
              "facets keep their programme_N labels and nothing is split off")
      derived <- setNames(rep(FALSE, ncol(sc)), names(sc))
    }

    sc$group <- factor(as.character(mouse_fib@meta.data[[grp]]), levels = lv)
    long <- sc %>% pivot_longer(-group, names_to = "programme",
                                values_to = "score")
    if (!is.null(nm))
      long$label <- factor(unname(nm[long$programme]), levels = unname(nm))
    else
      long$label <- factor(long$programme)
    long$derived <- unname(derived[long$programme])

    main <- long %>% filter(!derived)
    if (nrow(main))
      save_plot(box_panel(main, "group", "score", "label",
                          ylab = "Module score"),
                "explor_programme_scores", 6.6, 5.0)
    # As of 2026-09-29 basement_membrane is data_derived, so this branch is
    # live and that panel is drawn here, outside the manuscript figures.
    dd <- long %>% filter(derived)
    if (nrow(dd)) {
      message("  descriptive (data-derived) programmes, drawn separately: ",
              paste(unique(as.character(dd$label)), collapse = ", "))
      save_plot(box_panel(dd, "group", "score", "label",
                          ylab = "Module score"),
                "explor_programme_scores_data_derived", 4.4, 2.6)
    }
  } else {
    message("  no programme_* columns in the object; Fig3F skipped")
  }

  if (!is.null(CONFIG$human_caf_rds) && file.exists(CONFIG$human_caf_rds)) {
    ensure("human_caf")
    hd <- expr_long(human_caf, VLN_GENES_HUMAN) %>%
      mutate(group = factor(group, levels = levels(Idents(human_caf))))
    save_plot(violin_panel(hd, "group", "expression", "gene",
                           ylab = "Normalised expression") +
                theme(strip.text = element_text(face = "italic")),
              "explor_human_CAF_violin", 6.4, 5.6)
  }
  invisible(TRUE)
}


# =============================================================================
# ENRICHMENT DOT PLOTS
# -----------------------------------------------------------------------------
# clusterProfiler's dotplot() inherits the ambient theme and does not wrap term
# names, so long GO descriptions run into each other and the dots are sized for
# a full-page figure. Rebuilt here from the result table.
# =============================================================================

# "k/n" -> k and n as numbers.
ratio_parts <- function(s) {
  v <- strsplit(as.character(s), "/", fixed = TRUE)
  list(num = vapply(v, function(z) as.numeric(z[1]), numeric(1)),
       den = vapply(v, function(z) as.numeric(z[2]), numeric(1)))
}

# --- GO names that are too long to draw --------------------------------------
# ADDED 2026-10-07. At FIG$go_wrap = 24 one term in the Cxcl12-high panel,
# "antigen processing and presentation of endogenous peptide antigen via MHC
# class I", wrapped to FOUR lines where the row height allows three, and its
# last line physically overlapped the name of the term below it in the
# published PDF. The alternatives were fewer terms (FIG$go_n_terms 12 -> 10,
# which would have made the Results sentence "the twelve terms shown" wrong)
# or a taller panel; shortening one name costs neither.
#
# RULE: a replacement may only DROP words, never add or reinterpret them, and
# the term must still be identifiable. Table S3 carries every term under the
# name GO gave it, and the Figure 3D legend has to say so.
GO_DISPLAY_NAMES <- c(
  "antigen processing and presentation of endogenous peptide antigen via MHC class I" =
    "antigen presentation via MHC class I")

go_display <- function(x) {
  hit <- x %in% names(GO_DISPLAY_NAMES)
  if (any(hit)) {
    x[hit] <- GO_DISPLAY_NAMES[x[hit]]
    message("  shortened for display: ", sum(hit),
            " term name(s); Table S3 keeps the full name")
  }
  x
}

go_dotplot <- function(res, title = NULL, n = FIG$go_n_terms,
                       wrap = FIG$go_wrap, p_limits = NULL,
                       count_limits = NULL, log_colour = FALSE) {
  df <- as.data.frame(res)
  if (!nrow(df)) return(NULL)

  gr <- ratio_parts(df$GeneRatio)
  bg <- ratio_parts(df$BgRatio)
  df$ratio  <- gr$num / gr$den
  df$bg_n   <- bg$num                    # genes annotated to the term
  # clusterProfiler >= 4.12 supplies FoldEnrichment; older versions do not.
  if (!"FoldEnrichment" %in% names(df))
    df$FoldEnrichment <- df$ratio / (bg$num / bg$den)

  n_all <- nrow(df)
  keep <- df$Count >= FIG$go_min_count
  if (!is.na(FIG$go_max_bg)) keep <- keep & df$bg_n <= FIG$go_max_bg
  if (length(FIG$go_include)) keep <- keep | df$Description %in% FIG$go_include
  dropped <- df[!keep, , drop = FALSE]
  df <- df[keep, , drop = FALSE]
  if (!nrow(df)) {
    message("  ", title, ": no term passed the filters")
    return(NULL)
  }
  # Selection is always by adjusted P. Only the display order is configurable.
  df <- df %>% arrange(p.adjust) %>% head(n)
  if (length(FIG$go_include)) {
    forced <- setdiff(intersect(FIG$go_include, dropped$Description),
                      df$Description)
    if (length(forced))
      message("  forced in: ", paste(forced, collapse = "; "))
  }
  # What the filters cost, among terms that would otherwise have made the plot.
  near <- dropped %>% arrange(p.adjust) %>% head(5)
  if (nrow(near) && any(near$p.adjust < max(df$p.adjust)))
    message("  excluded by the filters but more significant than the weakest ",
            "term shown: ",
            paste(sprintf("%s (n=%d, bg=%d, p.adj %.1e)", near$Description,
                          near$Count, near$bg_n, near$p.adjust),
                  collapse = "; "))

  xvar <- if (identical(FIG$go_x, "fold")) "FoldEnrichment" else "ratio"
  xlab <- if (identical(FIG$go_x, "fold")) "Fold enrichment" else "Gene ratio"
  ord  <- if (identical(FIG$go_order, "x")) df[[xvar]] else -df$p.adjust

  message("  ", title, ": ", nrow(df), " of ", n_all,
          " terms shown (Count >= ", FIG$go_min_count,
          ", top ", n, " by adjusted P), ordered by ",
          if (identical(FIG$go_order, "x")) xlab else "adjusted P")
  chk <- df[order(-ord), c("Description", "Count", "bg_n", "FoldEnrichment",
                           "p.adjust")]
  for (i in seq_len(nrow(chk)))
    message(sprintf("    %-48s %3d/%4d  FE %5.2f  p.adj %8.1e",
                    substr(chk$Description[i], 1, 48), chk$Count[i],
                    chk$bg_n[i], chk$FoldEnrichment[i], chk$p.adjust[i]))

  # Shorten a handful of GO names FOR DISPLAY ONLY, before wrapping. Applied
  # here and nowhere else: selection, statistics and Table S3 all use the
  # Description as GO supplies it, so the shortening cannot change what is
  # shown, only how it is spelled on the axis. See GO_DISPLAY_NAMES.
  df$Description <- go_display(df$Description)
  df <- df %>%
    mutate(Description = str_wrap(Description, FIG$go_wrap),
           Description = factor(Description,
                                levels = Description[order(ord)]))

  ggplot(df, aes(.data[[xvar]], Description, size = Count,
                 colour = p.adjust)) +
    geom_point() +
    # p_limits and count_limits exist so that several panels can share one
    # guide: patchwork collects guides only when the scales are identical, and
    # without common limits each panel keeps its own legend. Three panels x two
    # legends is most of the width of a 180 mm figure.
    #
    # log_colour is what makes a shared P guide honest. Adjusted P here spans
    # 1e-18 to 1e-3, so on a linear scale the single most extreme term would
    # occupy the whole range and every other dot would read as the same colour.
    # On log10 the range is six even steps.
    scale_size_continuous(range = FIG$go_dot_range, name = "Genes",
                          limits = count_limits) +
    scale_colour_gradient(low = FIG$palette[1], high = "grey70",
                          name = "Adjusted P", limits = p_limits,
                          trans = if (log_colour) "log10" else "identity",
                          labels = function(v) format(v, digits = 1,
                                                      scientific = TRUE),
                          # ADDED 2026-10-07. The default bar is about 1 in
                          # wide and carries three 5-character labels, which
                          # need about 1.3 in: in the published PDF they
                          # overlapped into "1e-14e-10e-06". Widening the bar
                          # is the fix that keeps all three breaks and the
                          # scientific notation.
                          guide = guide_colourbar(
                            barwidth  = unit(1.5, "in"),
                            barheight = unit(0.10, "in"))) +
    labs(x = xlab, y = NULL, title = rich_title(title)) +
    theme_fig() +
    theme(axis.text.y = element_text(lineheight = 0.9),
          plot.title  = md_elem(size = FIG$title_size, hjust = 0),
          # ADDED 2026-10-07. The rightmost panel's title is longer than its
          # plotting area, so "Up in Meflin+/Thy1high" ran 4.9 pt past the
          # right edge of the page and lost its last word in the PDF.
          plot.margin = margin(5.5, 18, 5.5, 5.5))
}

# =============================================================================
# REMOVED 2026-09-23: the hypothesis-driven morphogen-antagonist panel
# -----------------------------------------------------------------------------
# MORPHOGEN_ANTAGONISTS and fig_antagonists() drew Figure S1F: secreted Wnt and
# hedgehog antagonists at gene level across the three subsets, with a legend
# calling them "pre-specified".
#
# They were not pre-specified. The list was assembled on 2026-09-14 from the top
# markers of the Thy1-high subset, so the panel showed that genes taken from one
# subset's marker list are restricted to that subset. Figure S1F was dropped on
# 2026-09-16 and the code is removed here so that the panel, and the S1D_04b /
# S1D_04c tables that went with it, cannot be regenerated.
#
# Notum, Wif1, Dkk2, Nkd1, Nkd2, Hhip and Kremen1 remain in the subset marker
# dot plot, where they are shown as markers -- which is what they are.
# =============================================================================

# =============================================================================
# DRAFT PANEL, ADDED 2026-10-05 -- Figure 3F candidate, not yet in fig_panels().
#
# One picture for the whole enrichment paragraph: the normalised enrichment
# score of three gene sets across the three subsets, one set standing for each
# subset's characterisation in the Results. Unlike Figure 3D it is
# threshold-free and signed, so "depleted" and "not enriched" are visible
# rather than absent.
#
# WHY THE SET LIST IS FIXED HERE AND NOT CHOSEN FROM THE OUTPUT. The
# one-versus-rest GSEA returns 221 significant sets. Picking a handful of them
# after seeing the result is the same mistake that cost Figure S1F: the panel
# would then be a drawing of the conclusion rather than evidence for it. The
# rule is therefore declared: ONE set per subset, and each one is the set
# NAMED IN THE RESULTS for that subset. Changing this list means changing the
# Results, and the Methods must state the rule if the panel is used.
#
# Non-significant values are drawn, hollow. That is the point of the panel --
# the Thy1-high subset's position on the matrix axis is a null, and a null
# cannot be shown by leaving the dot out.
# THE SET LIST AND THE RULE THAT PICKED IT. The one-versus-rest GSEA returns 221
# significant sets out of 12,733 tested. Choosing a handful after seeing them is
# the mistake that cost Figure S1F, so the rule is declared and the list follows
# from it:
#
#   (i)   the programme named in the Results for each subset -- one each;
#   (ii)  for the matrix axis, which carries the comparison, one corroborating
#         set from a SECOND collection, so the axis is not resting on one
#         annotation;
#   (iii) the one axis that separates the two Thy1-low subsets from each other;
#   (iv)  the confounder the Methods already name, because leaving it out of a
#         panel that shows everything else is a choice a reader cannot see.
#
# Everything else is in Table S4b. Adding a set here means changing the rule
# above and saying so in the Methods.
#
# ORDER, 2026-10-06: the rows follow the order in which the Results make the
# claims -- Cxcl12-high, then the matrix axis, then Thy1-high -- with the three
# matrix annotations adjacent and the confounder last. The first version listed
# them in the order the rule generated them, which put the two matrix sets at
# the top and the third (EMT) five rows below, so the reader had to find them.
GSEA_PANEL_SETS <- c(
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",            # (i)   Cxcl12-high
  "extracellular matrix organization",             # (i)   Mfap5-high
  # CHANGED 2026-10-06, was REACTOME_COLLAGEN_FORMATION. That set is the
  # NARROWER concept, and the choice of which narrower concept to show was
  # mine, not the rule's; it is also not significant in the Mfap5-high subset
  # (+1.6, open), so as a corroboration it corroborated only the depletion.
  # The Reactome set of the SAME NAME as the GO set is what rule (ii) asks for
  # -- one concept, two databases -- and it is significant in both directions
  # (Mfap5-high +1.88, Cxcl12-high -1.75).
  "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",    # (ii)  the same axis, Reactome
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",    # (iii) Mfap5-high vs Cxcl12-high
  "negative regulation of Wnt signaling pathway",  # (i)   Thy1-high
  "HALLMARK_MYC_TARGETS_V1")                       # (iv)  the depth confounder

# Axis labels. The Description strings are the databases' own: one shouted in
# upper case with a prefix, the others lower case. Printing them raw makes the
# panel look like a console dump. Selection is by the names above and display
# by these, so renaming a label can never change which set is drawn.
GSEA_PANEL_LABELS <- c(
  "extracellular matrix organization"            = "Extracellular matrix organisation (GO)",
  "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION"   = "Extracellular matrix organisation (Reactome)",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE"           = "Interferon-\u03b1 response (hallmark)",
  "negative regulation of Wnt signaling pathway" = "Wnt signalling, negative regulation (GO)",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"   = "Epithelial-mesenchymal transition (hallmark)",
  "HALLMARK_MYC_TARGETS_V1"                      = "MYC targets V1 (hallmark)")

# layout = "dodged": one axis, subsets as colour. The comparison the paragraph
#   makes is BETWEEN subsets on the same axis (+2.09 against -1.96), and that
#   comparison is read off one axis, not across three facets.
# layout = "facet": the earlier three-panel version, kept so the two can be
#   put side by side before one is chosen.
fig_named_gene_violins <- function(sets = NAMED_GENE_SETS,
                                   compact = NAMED_GENE_COMPACT,
                                   set_labels = NAMED_GENE_BLOCK_STRIP) {
  message("Gene-level panel (Figure S1D)",
          if (compact) paste0(" [compact",
                              if (!set_labels) ", gene names only" else "", "]")
          else " [full]")
  ensure("mouse_fib")
  grp <- if ("subset_display" %in% colnames(mouse_fib@meta.data))
    "subset_display" else "subset"
  lv <- levels(mouse_fib@meta.data[[grp]])
  genes <- unlist(sets, use.names = FALSE)

  d <- expr_long(mouse_fib, genes, grp)        # reports absent genes itself
  got <- unique(as.character(d$gene))
  if (!length(got)) return(invisible(NULL))
  if (length(setdiff(genes, got)))
    message("  NOT DRAWN (absent from the object): ",
            paste(setdiff(genes, got), collapse = ", "),
            " -- the Results name these genes, so check that sentence")

  # Facet titles carry which signature each gene belongs to, because the point
  # of the panel is the contrast between the two sets.
  which_set <- setNames(rep(names(sets), lengths(sets)), genes)
  d$gene <- factor(d$gene, levels = intersect(genes, got))
  d$group <- factor(d$group, levels = lv)
  # The gene name is italic and the block label is not. The strip therefore
  # has to be rendered as markdown; with ggtext absent, use_rich() is FALSE and
  # the whole strip falls back to plain text rather than printing the asterisks.
  rich <- use_rich()
  g <- levels(d$gene)
  # The block name, marked up when ggtext is there. NAMED_GENE_BLOCK_LABELS holds the
  # marked-up form for the blocks that need one; the rest are set as written.
  blk <- which_set[g]
  if (rich) {
    hit <- blk %in% names(NAMED_GENE_BLOCK_LABELS)
    blk[hit] <- NAMED_GENE_BLOCK_LABELS[blk[hit]]
  }
  lab <- if (!set_labels) {
    # What each gene is has to be named in the figure legend when it is not on
    # the strip -- see the note at NAMED_GENE_BLOCK_STRIP.
    setNames(if (rich) paste0("*", g, "*") else g, g)
  } else if (rich) {
    setNames(paste0("*", g, "*<br>(", blk, ")"), g)
  } else {
    setNames(paste0(g, "\n(", blk, ")"), g)
  }
  strip_sz <- if (!compact) FIG$strip_size
              else if (set_labels) FIG$strip_size - 1 else FIG$strip_size - 0.5

  p <- violin_panel(d, "group", "expression", "gene", ncol = NAMED_GENE_NCOL,
                    ylab = "Normalised expression",
                    x_axis = !compact, legend = compact) +
    facet_wrap(~ gene, ncol = NAMED_GENE_NCOL, scales = "free_y",
               labeller = labeller(gene = lab)) +
    theme(strip.text = if (!compact) {
            # unchanged from the first drawn version, so compact = FALSE
            # reproduces that panel exactly
            if (rich) md_elem(size = strip_sz, colour = "black",
                              lineheight = 1.1)
            else element_text(size = strip_sz, colour = "black")
          } else if (rich) {
            md_elem(size = strip_sz, colour = "black", lineheight = 1.0,
                    margin = margin(b = 1.2, t = 0.5))
          } else {
            element_text(size = strip_sz, colour = "black", lineheight = 0.95,
                         margin = margin(b = 1.2, t = 0.5))
          })
  if (!compact) return(p)

  # Three y breaks per facet, not the default five. With fifteen free y scales
  # the tick labels are a real share of the width, and the panel is read for
  # which subset is highest, not for the value.
  p +
    scale_y_continuous(n.breaks = 3, expand = expansion(mult = c(0.02, 0.06))) +
    theme(axis.text.y      = element_text(size = FIG$base_size - 1.5),
          axis.title.y     = element_text(size = FIG$base_size),
          legend.text      = md_elem(size = FIG$base_size - 1,
                                     colour = "black"),
          legend.key.size  = unit(3, "mm"),
          panel.spacing.x  = unit(0.06, "in"),
          panel.spacing.y  = unit(0.09, "in"),
          plot.margin      = margin(2, 3, 1, 1))
}

fig_gsea_nes <- function(sets = GSEA_PANEL_SETS, alpha = 0.05,
                         labels = GSEA_PANEL_LABELS,
                         layout = c("dodged", "facet")) {
  layout <- match.arg(layout)
  message("GSEA NES panel (draft Figure 3F, ", layout, ")")
  f <- tab("Fig5_12c_GSEA_one_vs_rest_all.csv")
  if (!file.exists(f)) {
    message("  ", basename(f), " not found. Run run_steps(14) once; the ",
            "GSEA step writes it and caches the result, so later runs reuse ",
            "the cache (CONFIG$redo_gsea <- TRUE forces a recompute). The ",
            "older Fig5_12_GSEA_one_vs_rest.csv holds significant sets only ",
            "and therefore cannot show a null.")
    return(invisible(NULL))
  }
  d <- utils::read.csv(f, stringsAsFactors = FALSE)
  d <- d[d$Description %in% sets,
         c("subset", "collection", "Description", "NES", "p.adjust")]
  missing <- setdiff(sets, d$Description)
  if (length(missing))
    message("  not found in the GSEA output (check the spelling against the ",
            "Description column): ", paste(missing, collapse = "; "))
  if (!nrow(d)) return(invisible(NULL))

  # A name can occur in more than one collection, and two dots on one row would
  # read as two measurements of the same thing. Refuse to guess: name the clash
  # and keep the first, so a silent double dot cannot reach a figure.
  dup <- duplicated(d[, c("subset", "Description")])
  if (any(dup)) {
    message("  the same set name appears in more than one collection (",
            paste(unique(d$Description[dup]), collapse = "; "),
            "); keeping the first.")
    d <- d[!dup, , drop = FALSE]
  }

  # Subset order comes from the object, as every other panel's does. An earlier
  # version tested for a SUBSET constant that does not exist in this script and
  # fell through to sort(), which put the Thy1-high subset second.
  lv <- tryCatch({ ensure("mouse_fib")
                   levels(idents_of(mouse_fib, "subset")) },
                 error = function(e) NULL)
  lv <- if (is.null(lv)) sort(unique(d$subset)) else
    intersect(lv, unique(d$subset))

  grid <- expand.grid(subset = lv, Description = sets, stringsAsFactors = FALSE)
  d <- merge(grid, d, by = c("subset", "Description"), all.x = TRUE)
  n_gap <- sum(is.na(d$NES))
  if (n_gap)
    message("  ", n_gap, " subset x set combination(s) absent from the GSEA ",
            "output; drawn as gaps, never as zeros")

  disp <- display_of(lv)
  lab  <- ifelse(sets %in% names(labels), labels[sets], sets)
  d$label <- factor(d$Description, levels = rev(sets),
                    labels = rev(str_wrap(lab, 26)))
  d$sig   <- !is.na(d$p.adjust) & d$p.adjust < alpha
  d       <- d[!is.na(d$NES), , drop = FALSE]

  if (layout == "facet") {
    d$subset <- factor(d$subset, levels = lv, labels = str_wrap(disp, 14))
    p <- ggplot(d, aes(x = NES, y = label)) +
      geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey70") +
      geom_point(data = d[d$sig, , drop = FALSE], aes(fill = NES),
                 shape = 21, size = 2.6, stroke = 0.4, colour = "grey25") +
      geom_point(data = d[!d$sig, , drop = FALSE], shape = 21, size = 2.6,
                 stroke = 0.4, fill = NA, colour = "grey55") +
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                           midpoint = 0, name = "NES") +
      facet_wrap(~ subset, nrow = 1) +
      # 6 pt was not enough: the Cxcl12-high panel ends at +2.45 and the
      # Mfap5-high panel begins at -2.42, so two unrelated dots sat side by
      # side across the boundary and read as a pair.
      theme(panel.spacing.x = unit(14, "pt"))
  } else {
    # CHANGED 2026-10-07: was scales::hue_pal(), chosen to match what DimPlot
    # happened to draw in Figure 3A. umap_plot() now passes FIG$palette to
    # DimPlot instead, so matching 3A means using the same list -- and the
    # three panels that carry subset colour (3A2, 3E, 3F) finally agree.
    cols <- setNames(FIG$palette[seq_along(lv)], disp)
    d$subset <- factor(d$subset, levels = lv, labels = disp)
    # reverse = TRUE because the y scale is reversed: without it the three
    # subsets run bottom-to-top within each row while the key runs top-to-
    # bottom, so the key and the rows disagree about the order. The argument
    # arrived in ggplot2 3.5.0; on anything older the panel is drawn as before
    # rather than failing, and the message says which it was.
    pos <- if ("reverse" %in% names(formals(position_dodge)))
      position_dodge(width = 0.62, reverse = TRUE)
    else {
      message("  ggplot2 < 3.5: the within-row order stays the reverse of the ",
              "key (upgrade ggplot2 to fix)")
      position_dodge(width = 0.62)
    }
    p <- ggplot(d, aes(x = NES, y = label, colour = subset, group = subset)) +
      geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey70") +
      # The segment is what makes the direction readable at a glance; the dot
      # alone leaves the reader measuring against the zero line.
      geom_linerange(aes(xmin = 0, xmax = NES), position = pos,
                     linewidth = 0.5) +
      geom_point(aes(fill = ifelse(sig, as.character(subset), NA_character_)),
                 shape = 21, size = 2.4, stroke = 0.6, position = pos) +
      # labels = lab_fn() with a markdown legend.text: without it this was the
      # only panel in Figures 3 and S1 printing the subset names as stored,
      # "Meflin+/Cxcl12-high", while every other panel set them with real
      # superscripts.
      scale_colour_manual(values = cols, name = NULL, labels = lab_fn()) +
      scale_fill_manual(values = cols, na.value = "white", guide = "none") +
      guides(colour = guide_legend(reverse = FALSE))
  }

  p +
    # n.breaks 5 -> 3: on a -3 to 3 range five breaks gave seven labels and
    # "-3-2-1" ran together into one string at this panel width.
    scale_x_continuous(n.breaks = 3,
                       expand = expansion(mult = c(0.12, 0.12))) +
    labs(x = "Normalised enrichment score", y = NULL,
         caption = paste0("filled, adjusted P < ", alpha,
                          "; open, tested and not significant")) +
    theme_fig() +
    # legend.text LAST. theme_fig() sets it with element_text(), so setting it
    # inside the dodged branch above (which is added BEFORE theme_fig()) was
    # silently overwritten and the 2026-10-07 panel printed the key as
    # "Meflin<sup>+</sup>/Cxcl12<sup>high</sup>". A later theme() wins, so the
    # markdown element has to come after theme_fig(), not before it.
    theme(plot.caption = element_text(hjust = 0, size = rel(0.85)),
          legend.text  = md_elem(size = FIG$base_size, colour = "black"))
}

fig_go <- function() {
  message("Enrichment dot plots")
  f <- obj("05_go_results.rds")
  if (!file.exists(f)) {
    message("  ", basename(f), " not found; run section 13 first")
    return(invisible(FALSE))
  }
  go <- readRDS(f)

  # CHANGED 2026-09-23, again 2026-09-24. Figure S1D was the pairwise
  # peribiliary-versus-periportal enrichment; it became the one-versus-rest
  # enrichment, and that panel is now Figure 5D, in the main figure. The
  # enrichment is the conclusion of the single-cell analysis -- the sentence the
  # Results build to is that the subsets have distinguishable programmes -- so it
  # belongs with the panels that define the subsets, not in the supplement.
  # Figure S1D is now the marker violin plot, moved out of Figure 5.
  # The pairwise contrast described two subsets out of three, less strongly
  # (best adjusted P 1.2e-6 against 8.9e-9; ECM organization fold enrichment 6.50
  # against 7.14), and never mentioned the Cxcl12-high subset, which is 65 % of
  # the fibroblasts. Nothing the pairwise panel showed is absent from this one.
  # The pairwise DEG table stays in the supplement; only the figure was dropped.
  #
  # Each panel keeps its own colour legend. The adjusted P ranges differ by
  # orders of magnitude between subsets, so a collected guide would be dominated
  # by the single most extreme term and every other dot would read as the same
  # colour. Two legends cost width; one shared legend would misrepresent a panel.
  if (!is.null(go$subsets)) {
    # Common scales across the three panels: patchwork collects guides only
    # when the scales match, and three panels x two guides is most of the
    # width. The P scale is log10 because the range spans 1e-18 to 1e-3.
    shown <- Filter(Negate(is.null), lapply(go$subsets, function(r) {
      if (is.null(r)) return(NULL)
      d <- as.data.frame(r$simplified)
      d <- d[d$Count >= FIG$go_min_count, , drop = FALSE]
      if (!nrow(d)) return(NULL)
      d[order(d$p.adjust), ][seq_len(min(FIG$go_n_terms, nrow(d))), ]
    }))
    plim <- range(unlist(lapply(shown, `[[`, "p.adjust")), na.rm = TRUE)
    nlim <- range(unlist(lapply(shown, `[[`, "Count")), na.rm = TRUE)
    sp <- lapply(names(go$subsets), function(cl) {
      r <- go$subsets[[cl]]
      # CHANGED 2026-10-07: was paste("Up in", display_of(cl)). The title is
      # left-aligned at the panel's left edge, so a title longer than the panel
      # overhangs to the right no matter how much right margin is reserved --
      # which is why adding plot.margin did not fix the third panel losing
      # "1high" off the page. Dropping "Up in " is 22 pt, four times the 4.9 pt
      # overhang, and it makes the titles read exactly like the subset labels
      # in 3B, 3E, 3F and S1D. THE DIRECTION NOW HAS TO BE IN THE LEGEND:
      # "terms over-represented among the genes upregulated in each subset".
      if (is.null(r)) NULL else
        go_dotplot(r$simplified, display_of(cl),
                   p_limits = plim, count_limits = nlim, log_colour = TRUE)
    })
    sp <- Filter(Negate(is.null), sp)
    if (length(sp)) {
      w <- sz("go")
      if (identical(FIG$go_layout, "column")) {
        # Stacked, so the x axis is drawn once and the guides collect on the
        # right, where a vertical stack has room.
        for (i in seq_len(length(sp) - 1L))
          sp[[i]] <- sp[[i]] +
            theme(axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.title.x = element_blank())
        save_plot(wrap_plots(sp, ncol = 1, guides = "collect") &
                    scale_x_continuous(
                      expand = expansion(mult = c(0.18, 0.12))) &
                    theme(legend.position = FIG$go_legend,
                          legend.box = if (identical(FIG$go_legend, "bottom"))
                            "horizontal" else "vertical",
                          legend.key.height = unit(2.5, "mm"),
                          legend.key.width  = unit(4, "mm")),
                  "Fig3D_GO", w[1], w[2])
      } else {
        # Side by side. Each panel keeps its own y axis, because the three
        # subsets share almost no terms, and its own x axis, because the scales
        # are free. The guides are collected once underneath.
        #
        # n.breaks = 3 is here because the plotting area is under an inch wide
        # (see PANEL_SIZE$go): the default four or five breaks put the tick
        # labels close enough to touch at base_size 8.
        #
        # The expansion is the fix for the dot that sat on the y axis. ggplot's
        # default for a continuous scale is 5 % of the range at each end, which
        # is less than the radius of a size-4 dot in a panel this narrow, so
        # the term with the lowest fold enrichment -- always present, because
        # the panel is sorted by it -- was drawn half outside the panel and
        # clipped. 18 % on the left clears the largest dot; 12 % on the right
        # keeps the top term off the border without wasting width.
        save_plot(wrap_plots(sp, nrow = 1, guides = "collect") &
                    scale_x_continuous(
                      n.breaks = 3,
                      expand = expansion(mult = c(0.18, 0.12))) &
                    theme(legend.position = "bottom",
                          legend.box = "horizontal",
                          legend.key.height = unit(2.5, "mm"),
                          legend.key.width  = unit(4, "mm"),
                          panel.spacing.x = unit(2, "mm"),
                          plot.margin = margin(2, 3, 2, 2)),
                  "Fig3D_GO", w[1], w[2])
      }
    }
  } else {
    message("  no per-subset GO in 05_go_results.rds; re-run section 13")
  }

  # REMOVED 2026-09-23: FigS1E_GSEA_GO. Across GO biological process and three
  # pre-specified MSigDB collections the only significant gene set was
  # HALLMARK_MYC_TARGETS_V1, whose leading edge is initiation and elongation
  # factors, chaperonins, spliceosome and proteasome components -- the
  # translation machinery. The subsets differ in ribosomal transcript content
  # (17.0 / 13.5 / 12.8 %), and the ^Rp[sl] exclusion removes the ribosomal
  # proteins only, so this is the depth confounder under a biological heading.
  # The GO BP result was 20 terms at adjusted P 0.013-0.044, three of them the
  # same five genes counted as eye terms. The analysis is reported as a negative
  # result in the Methods and drawn nowhere.
  invisible(TRUE)
}


# =============================================================================
# DRIVERS
# =============================================================================

fig_mouse <- function() { fig_umaps(); fig_dotplots(); fig_violins(); fig_go() }
fig_human <- function() { fig_umaps(); fig_dotplots(); fig_violins() }

# fig_all() writes the manuscript panels FIRST. It did not, for one version,
# and that is why three requested changes to the feature plots (then Figure 5B,
# now Figure 5C) appeared to do nothing: only fig_panels() writes them,
# fig_all() never called it, and a stale file
# from an older run sat in the folder looking current. The exploratory
# functions no longer write any manuscript filename either -- two writers for
# one file, with different settings, is the same defect that once labelled the
# panels with the internal codes PBF / PPF / CxF.
fig_all <- function(exploratory = TRUE) {
  fig_panels()
  if (exploratory) {
    message("\nExploratory panels (explor_*, not for the manuscript)")
    fig_umaps(); fig_dotplots(); fig_violins()
  }
  invisible(TRUE)
}

message("figures_final.R ", FIG_VERSION, " loaded from ",
        tryCatch(normalizePath("figures.R"), error = function(e) "?"))
message("fig_panels() writes every manuscript panel at its final printed size.",
        "\n  Exploratory: fig_umaps() / fig_dotplots() / fig_violins() / fig_go().  Also name_check() / caf_check() / caf_markers() / coexpr_check() / sample_check() / resolution_check() / trace_cells().")
message("Settings are in the FIG list; change one and re-run that function.")


# =============================================================================
# PANEL SIZES
# -----------------------------------------------------------------------------
# The composites are gone. fig5() and figS1() used to assemble the panels into
# one page; the figures are laid out in Illustrator instead, so all this script
# has to do is write each panel at the size it will occupy there.
#
# That removes the entire class of problem the composites created. Every awkward
# thing in the last week -- panels letterboxing at a fixed aspect, the column
# arithmetic, the enrichment panels squeezed until the dots vanished, the
# feature plots at a quarter size -- existed only because panels shared a page.
# Two of them were outright bugs: a second set of point sizes that drifted out
# of step with the first, and a feature_plot call missing group_by, which
# labelled the panels with the internal codes PBF / PPF / CxF. Neither could
# have happened with one panel per file.
#
# ONE RULE: each panel is written at its final printed size, so Illustrator
# never has to scale it. Scaling is not free here -- the scatter layers are
# rasterised, so enlarging degrades them, and scaling a panel changes the
# effective size of every font on it. Place at 100 %.
#
# The sizes below assume a 180 mm (7.09 in) full-width figure. If a panel ends
# up somewhere else, change its entry in PANEL_SIZE and re-run; the dot radius
# follows the width automatically through pt_for(..., width_in).
PANEL_SIZE <- list(
  # Figure 5
  fib_umap        = c(3.4, 2.6),   # numbered key adds a legend column
  allcell_umap_hl = c(2.6, 2.6),   # highlighted, no legend
  # CHANGED 2026-09-29: was c(4.8, 5.3) for 4 genes in 2 x 2. Six genes in
  # 3 x 2, sub-plot width held at 2.4 in so the dot radius, and therefore
  # Fig3A2, do not have to be retuned.
  fib_features    = c(7.2, 5.3),   # 3 x 2, each sub-plot about 2.4 in square
  # CHANGED 2026-09-29: c(7.0, 3.2) -> c(4.6, 4.8) -> c(5.0, 2.7).
  # The 7-facet strip lost myofibroblast, apCAF and iCAF (this is not a tumour
  # and there are no CAFs in it to classify), which left four and a 2 x 2. It
  # then lost basement_membrane to the data_derived flag, which left three:
  # universal_Pi16, universal_Col15a1 and adventitial, in one row.
  #
  # CHANGED 2026-09-29: adventitial is now data_derived too (its four genes are
  # markers of the subset it scores highest in, and no sentence in the Results
  # used it), so Fig. 3E is the two externally defined universal programmes.
  # c(5.0, 2.7) -> c(3.6, 2.7): the facets keep their width, one is removed.
  programme       = c(3.6, 2.7),   # 1 x 2 square facets + angled labels
  fib_dotplot     = c(2.6, 5.0),   # flipped: 24 genes down the side
  # ADDED 2026-10-06, Figure 3E since 2026-10-08: six gene sets down the side,
  # three subsets
  # dodged within each row. Width is set by the longest set name after
  # wrapping at 26 characters, not by the data.
  gsea_nes        = c(4.8, 2.6),
  # Figure S1D: 18 genes, 6 columns, 3 rows. Three sizes, picked by
  # NAMED_GENE_COMPACT / NAMED_GENE_BLOCK_STRIP; named_gene_size() below does the picking.
  #
  # CHANGED 2026-10-07: c(6.8, 5.0) -> c(5.8, 3.55). The first drawn version
  # did not fit on the Figure S1 page. The 1.45 in of height that went is
  # almost all the angled subset labels repeated under fifteen facets, now one
  # key under the panel; the width follows, because the facets no longer have
  # to be as wide as the longest angled label. Nothing was made smaller except
  # the strip (8 -> 7 pt) and the y tick labels (8 -> 6.5 pt), both still above
  # the 5 pt floor Cancer Research sets for final-size text.
  # CHANGED 2026-10-08 with the interferon block: 14 genes in 5 columns ->
  # 18 in 6 (NAMED_GENE_NCOL). One more column at the ~1.1 in each the five already
  # had: 5.8 -> 6.9 in, which is 175 mm, just inside the 180 mm full width.
  # The height does not move -- it is still three rows.
  named_gene_vln       = c(6.9, 3.55),  # compact, strip names what the gene is
  named_gene_vln_nolab = c(5.6, 3.15),  # compact, gene names only
  named_gene_vln_full  = c(8.0, 5.00),  # as first drawn, labels under facets
  # Two thirds of the previous width. Narrowing a GO dot plot is not free: the
  # panel is term text + plotting area + guides, so taking an inch and a half
  # off the width has to come out of one of the three. The text gives up the
  # least -- the wrap drops from 30 to 24 characters -- and the guides move
  # underneath, where they cost height instead of width. What is left for the
  # dots is about 2.3 in, more than the 1.4 in they would get with the guides
  # still on the right at this width.
  #
  # The height rises because narrower wrapping means more lines per term. Set
  # FIG$go_n_terms to 5 to get back to about 5.4 in.
  #
  # CHANGED 2026-09-29: was c(3.2, 6.2) stacked, then c(7.1, 6.0) side by side,
  # now c(7.8, 4.5) -- three quarters of that height and 10 % more width.
  # The three subsets are side by side (FIG$go_layout), so the height carries
  # one subset's 12 terms instead of all 36.
  #
  # The height is the tight dimension now. Of the 4.5 in, roughly 0.25 goes to
  # the panel title, 0.40 to the x axis, 0.45 to the collected guides and 0.15
  # to margins, leaving about 3.25 in for twelve terms, or 0.27 in each. A
  # two-line name at base_size 7 and lineheight 0.9 is about 0.18 in, so it
  # fits; a three-line name is about 0.26 in and only just fits. That is why
  # FIG$go_wrap went to 32. IF THE NAMES TOUCH AGAIN, the first thing to change
  # is FIG$go_n_terms (12 -> 10 buys 0.05 in per term), not the height.
  #
  # 7.8 in is 198 mm, wider than the 180 mm full-width figure. The panel is
  # assembled by hand, so it will be scaled there; at 180 mm everything above
  # shrinks by 9 %, which the term text can take and the dots do not care about.
  # CHANGED 2026-09-30: c(7.8, 4.5) -> c(6.6, 5.0). 15 % narrower, and taller
  # because 8 pt text wrapped at 24 characters puts more terms on three lines.
  # Rough budget per panel: term column about 1.3 in, plotting area about
  # 0.7 in, ticks and spacing about 0.2 in.
  #
  # This is close to the limit of what three panels side by side can do. If it
  # still reads crowded the next lever is FIG$go_n_terms (12 -> 10), NOT the
  # width: twelve terms in 0.7 in of plotting area is already tight. Changing
  # it means changing "twelve" in the Methods and the Figure 3D legend.
  go              = c(6.6, 5.0),   # 3 subsets side by side, shared guides
  # Figure S1
  allcell_umap    = c(6.0, 3.8),   # labelled, 16-entry key
  # CHANGED 2026-10-08: c(7.5, 3.5) -> c(7.5, 3.9). The guides moved from the
  # right of the panel to a band above it (dot_plot legend = "top"), which
  # hands roughly 90 pt of width back to the 44 gene columns -- they go from
  # about 8.9 to 11 pt apart. The band costs about 0.35 in of height, and at
  # 3.5 in the sixteen rows closed up until the dots touched between rows, so
  # the height buys exactly what the band took. The WIDTH is deliberately not
  # increased: 7.5 in is already past the 180 mm (7.09 in) full-column width,
  # and the point of the change is to use the width the panel already has.
  allcell_dotplot = c(7.5, 3.9),   # 44 genes x 16 populations, guides on top
  allcell_features= c(6.6, 4.4),
  violin          = c(6.4, 5.0),   # 12 genes, 4 x 3
  fib_features_s1 = c(6.6, 4.4)    # 6 genes, the reference version
)
sz <- function(name) {
  v <- PANEL_SIZE[[name]]
  if (is.null(v)) stop("no size for panel '", name, "' in PANEL_SIZE")
  v
}

# Figure S1D has three sizes because it has three layouts. Keeping the choice
# here rather than at the save site means the switch at NAMED_GENE_COMPACT is the only
# thing to change: the size cannot be left behind on the layout it was tuned
# for, which is how a panel ends up with the right facets and the wrong page.
named_gene_size <- function(compact = NAMED_GENE_COMPACT, set_labels = NAMED_GENE_BLOCK_STRIP) {
  if (!compact) sz("named_gene_vln_full")
  else if (set_labels) sz("named_gene_vln")
  else sz("named_gene_vln_nolab")
}

# Programme-score table in long form, with the data-derived programmes flagged.
# Factored out of fig_violins() so that the standalone panel and the composite
# cannot drift apart.
# --- programme names as they are printed ---------------------------------------
# ADDED 2026-10-07. Fig5_07_programme_definitions.csv holds the analysis
# script's internal names, and Figure 3E was drawing them straight onto the
# facet strips: "universal_Pi16" and "universal_Col15a1", underscore and all,
# for the two things the Results and Figure S1D call "Pi16+ universal" and
# "Col15a1+ universal". Cosmetic only, like GO_DISPLAY_NAMES: a name with no
# entry here is left exactly as the CSV spells it.
PROGRAMME_DISPLAY_NAMES <- c(
  universal_Pi16     = "Pi16+ universal",
  universal_Col15a1  = "Col15a1+ universal")

programme_display <- function(x) {
  x <- as.character(x)
  hit <- x %in% names(PROGRAMME_DISPLAY_NAMES)
  x[hit] <- PROGRAMME_DISPLAY_NAMES[x[hit]]
  x
}

programme_long <- function(x = mouse_fib, group_col = NULL) {
  grp <- if (!is.null(group_col)) group_col else
    if ("subset_display" %in% colnames(x@meta.data)) "subset_display" else "subset"
  lv <- levels(x@meta.data[[grp]])
  sc <- x@meta.data %>% select(matches("^programme_[0-9]+$")) %>% as_tibble()
  if (!ncol(sc)) return(NULL)
  defs <- tab("Fig5_07_programme_definitions.csv")
  if (file.exists(defs)) {
    dd <- read.csv(defs, stringsAsFactors = FALSE)
    nm      <- setNames(dd$programme, paste0("programme_", seq_len(nrow(dd))))
    derived <- setNames(as.logical(dd$data_derived),
                        paste0("programme_", seq_len(nrow(dd))))
  } else {
    message("  Fig5_07_programme_definitions.csv not found; ",
            "facets keep their programme_N labels and nothing is split off")
    nm <- NULL
    derived <- setNames(rep(FALSE, ncol(sc)), names(sc))
  }
  sc$group <- factor(as.character(x@meta.data[[grp]]), levels = lv)
  long <- sc %>% pivot_longer(-group, names_to = "programme",
                              values_to = "score")
  # Display names, applied here and nowhere else. The CSV carries the analysis
  # script's internal names, and Figure 3E was publishing them verbatim:
  # "universal_Pi16" and "universal_Col15a1", complete with the underscore,
  # while the Results and Figure S1D call the same two things "Pi16+ universal"
  # and "Col15a1+ universal". Same rule as GO_DISPLAY_NAMES -- the mapping is
  # cosmetic, nothing upstream of it changes, and a name with no entry is left
  # exactly as the CSV spells it.
  # nm[] <- keeps the names and the order, so the lookup below is unchanged.
  if (!is.null(nm)) nm[] <- vapply(nm, programme_display, character(1))
  long$label <- if (!is.null(nm))
    factor(unname(nm[long$programme]), levels = unname(nm)) else
    factor(long$programme)
  long$derived <- unname(derived[long$programme])
  long
}

# The whole-tissue UMAP as it appears in panel A. Nothing is removed from the
# object -- the grey cells are drawn, so the reader sees the extracted
# population as part of the tissue rather than as a plot of its own.
umap_highlight <- function(x, keep, mode = "single", group_by = "celltype",
                           grey = FIG$highlight_grey, width_in = NULL,
                           pt = NULL, label = FIG$highlight_label) {
  ids  <- droplevels(idents_of(x, group_by))
  miss <- setdiff(keep, levels(ids))
  if (length(miss))
    warning("umap_highlight: not a cluster in this object: ",
            paste(miss, collapse = ", "),
            ". Check FIG$highlight_keep against levels(mouse_all$celltype).",
            call. = FALSE)
  keep <- intersect(keep, levels(ids))
  if (!length(keep)) stop("umap_highlight: none of FIG$highlight_keep exists")

  sel <- as.character(ids) %in% keep
  message("  highlighting ", sum(sel), " of ", length(sel), " cells (",
          paste(keep, collapse = ", "), "), mode: ", mode)

  if (identical(mode, "single")) {
    lab  <- "extracted"
    col1 <- FIG$highlight_colour %||% FIG$palette[1]
    lv   <- c("other", lab)
    x$.hl <- factor(ifelse(sel, lab, "other"), levels = lv)
    cols  <- setNames(c(grey, col1), lv)
    shown <- lab
  } else {
    lv <- c("other", keep)
    x$.hl <- factor(ifelse(sel, as.character(ids), "other"), levels = lv)
    cols  <- setNames(c(grey, FIG$palette[seq_along(keep)]), lv)
    shown <- keep
  }
  if (is.null(pt)) pt <- pt_for(ncol(x), width_in = width_in)
  p <- DimPlot(x, reduction = umap_of(x), group.by = ".hl", label = FALSE,
               pt.size = pt, raster = FIG$raster,
               raster.dpi = FIG$raster_px, order = shown) +
    scale_colour_manual(values = cols, breaks = shown) +
    theme_fig() + umap_axes() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          legend.position = "none") +
    ggtitle(NULL) + aspect_layer()

  # Numerals from the FULL cluster set, so they agree with Fig. S1A, which is
  # where the reader finds out what each number is. "keep" draws only the
  # extracted clusters, which is legible on a small panel but leaves the grey
  # populations anonymous; "number" labels all of them.
  if (!identical(label, "none")) {
    lab <- cluster_label_layer(
      x, mode = "number", ids = ids, quiet = TRUE,
      size = FIG$feature_label_size,
      only = if (identical(label, "keep")) keep else NULL)
    if (!is.null(lab)) {
      ll <- if (inherits(lab, "gg") || inherits(lab, "Layer")) list(lab) else lab
      for (l in ll) p <- p + l
    }
  }
  p
}


# =============================================================================
# DRIVER FOR THE MANUSCRIPT PANELS
# =============================================================================

# Everything that goes into Figure 5 and Figure S1, each at its final size.
# --- one panel at a time ------------------------------------------------------
# ADDED 2026-10-08. fig_panels() redraws every manuscript panel, several of
# them feature plots over thousands of cells, which is a long wait when the
# thing being iterated on is one dot plot. Each panel is now its own function,
# and fig_panels() CALLS those functions rather than keeping its own copy of
# the code, so a panel drawn on its own and the same panel drawn by
# fig_panels() cannot drift apart. The loop while tuning a layout is
#
#     source("figures_final.R"); panel_s1c()
#
# or, for several, panels("3B", "S1D"). list_panels() prints the names.
#
# ensure() loads the .rds into the global environment on first use and reuses
# it after that, so only the FIRST call in a session pays for the read.
# save_plot() writes one file per panel: nothing else in figures_final/ is
# touched.
#
# WHEN ADDING A PANEL: write panel_xxx(), add it to PANEL_INDEX, and make
# fig_panels() call it. Never leave a second copy of the drawing code behind.

# The identity column the fibroblast panels group by. Computed here rather than
# once at the top of fig_panels() so that a panel drawn on its own groups the
# same way as one drawn in the full run.
fib_grp <- function() {
  ensure("mouse_fib")
  if ("subset_display" %in% colnames(mouse_fib@meta.data))
    "subset_display" else "subset"
}

# Figure 3C's dot geometry, which 3A2 and 3A1 are sized to match: the
# fibroblast UMAP and the feature plots sit side by side and hold the same
# cells, so their dots have to PRINT at the same size, and they are different
# panel widths, so the canvas radius has to differ to achieve that.
#
# This used to live in fig_panels()'s local variables, which made the match a
# property of the CALL ORDER -- 3C had to be drawn before 3A2 or 3A2 got
# whatever was in scope. It is a pure function of PANEL_SIZE, FIG and the cell
# count, so computing it on demand costs nothing and cannot go stale.
fib_dot_target <- function() {
  ensure("mouse_fib")
  wB   <- sz("fib_features")
  subB <- wB[1] / FIG$fib_feature_ncol
  ptB  <- pt_for(ncol(mouse_fib), FIG$feature_point, FIG$feature_point_range,
                 width_in = subB, ref_in = FIG$feature_point_ref_in,
                 what = "Fig3C")
  list(wB = wB, subB = subB, ptB = ptB, tgt = dot_diameter_in(ptB, subB))
}

# --- Figure 3 -----------------------------------------------------------------

panel_3c <- function() {
  g <- fib_dot_target()
  save_plot(feature_plot(mouse_fib, FIG$fib_feature_genes, group_by = fib_grp(),
                         ncol = FIG$fib_feature_ncol, width_in = g$wB[1]),
            "Fig3C_featureplots_fibroblasts", g$wB[1], g$wB[2])
}

panel_3a2 <- function() {
  g   <- fib_dot_target()
  w   <- sz("fib_umap")
  ptA <- pt_for_diameter(g$tgt, w[1])
  ptA <- min(max(ptA, FIG$umap_point_range[1]), FIG$umap_point_range[2])
  message(sprintf("  Fig3A2 matched to Fig3C: printed dot %.4f in; ", g$tgt),
          sprintf("radius %.1f px at %.2f in (5B: %.1f px at %.2f in)",
                  ptA, w[1], g$ptB, g$subB))
  save_plot(umap_plot(mouse_fib, group_by = fib_grp(), key = FIG$fib_umap_key,
                      width_in = w[1], pt = round(ptA, 1)),
            "Fig3A2_umap_fibroblasts", w[1], w[2])
}

# Sized from the same reference as 3A2, at a fraction of it.
panel_3a1 <- function() {
  ensure("mouse_all")
  tgt <- fib_dot_target()$tgt
  w1  <- sz("allcell_umap_hl")
  pt1 <- pt_for_diameter(tgt * FIG$highlight_dot_frac, w1[1])
  pt1 <- min(max(pt1, FIG$umap_point_range[1]), FIG$umap_point_range[2])
  message(sprintf("  Fig3A1: printed dot %.4f in (%.0f%% of Fig3C), radius %.1f px",
                  dot_diameter_in(pt1, w1[1]), 100 * FIG$highlight_dot_frac, pt1))
  save_plot(umap_highlight(mouse_all, FIG$highlight_keep, width_in = w1[1],
                           pt = round(pt1, 1)),
            "Fig3A1_umap_allcells_highlight", w1[1], w1[2])
}

panel_3b <- function() {
  x <- mouse_fib; Idents(x) <- idents_of(x, fib_grp())
  w <- sz("fib_dotplot")
  save_plot(dot_plot(x, FIB_DOT_GENES, dot_scale = 3, flip = TRUE),
            "Fig3B_dotplot_fibroblast_subsets", w[1], w[2])
}

# fig_go() draws this one; the wrapper exists so that every manuscript panel
# can be named the same way.
panel_3d <- function() fig_go()

panel_3f <- function() {
  long <- programme_long()
  if (is.null(long)) return(invisible(NULL))
  main <- long %>% filter(!derived)
  if (nrow(main)) {
    w <- sz("programme")
    # The two externally defined universal programmes in one row, square
    # facets. PANEL_SIZE$programme is sized for two; say so if that changes
    # rather than silently drawing a row with a hole in it.
    n_prog <- nlevels(droplevels(factor(main$label)))
    if (n_prog != 2)
      message("  Fig3E: ", n_prog, " programmes, not 2; the 1 x 2 layout ",
              "and PANEL_SIZE$programme were set for two")
    save_plot(box_panel(main, "group", "score", "label",
                        ylab = "Module score", ncol = n_prog, aspect = 1),
              "Fig3F_programme_scores", w[1], w[2])
  }
  # Data-derived programmes are deliberately NOT written to a manuscript file
  # here. fig_violins() draws them as explor_programme_scores_data_derived.
  # Say which ones, so that a panel dropped from Fig. 5E is a decision on the
  # record rather than something that quietly stopped appearing.
  dd <- unique(as.character(long$label[long$derived]))
  if (length(dd))
    message("  Fig3E: not in the main panel (data-derived): ",
            paste(dd, collapse = ", "),
            " -- see explor_programme_scores_data_derived")
  invisible(TRUE)
}

# --- Figure S1 ----------------------------------------------------------------

panel_s1a <- function() {
  ensure("mouse_all")
  w <- sz("allcell_umap")
  save_plot(umap_plot(mouse_all, group_by = "celltype", width_in = w[1]),
            "FigS1A_umap_allcells", w[1], w[2])
}

panel_s1b <- function() {
  ensure("mouse_all")
  w <- sz("allcell_features")
  save_plot(feature_plot(mouse_all, FIG$allcell_feature_genes,
                         group_by = "celltype", ncol = 3, width_in = w[1]),
            "FigS1B_featureplots_allcells", w[1], w[2])
}

panel_s1c <- function() {
  ensure("mouse_all")
  y <- mouse_all; Idents(y) <- idents_of(y, "celltype")
  w <- sz("allcell_dotplot")
  save_plot(dot_plot(y, MOUSE_DOT_GENES, legend = "top"),
            "FigS1C_dotplot_allcells", w[1], w[2])
}

# ACCEPTED 2026-10-06 as Figure S1D, not Figure 3F -- Figure 3 already carried
# five panels and the single-cell section is the part of this manuscript a
# reviewer is most likely to call tangential, so the supplement looked like the
# place for a threshold-free corroboration.
#
# REVERSED 2026-10-08, ON LAYOUT: this is Figure 3E and the module scores moved
# to 3F. Figure S1 had no room left and Figure 3 had, and the new order is the
# order the Results read -- 3D over-representation, 3E the signed version of
# the same comparisons, 3F the published signatures -- where before the reader
# had to leave for the supplement in the middle of that paragraph and come
# back. The 2026-10-06 worry stands and is now a size worry: Figure 3 carries
# A1, A2, B, C, D, E and F. What it adds that Figure 3D cannot:
# over-representation can only show that a programme IS enriched, so the
# Results sentences about depletion in the Cxcl12-high subset and about the
# null in the Thy1-high subset had no figure at all before this one.
panel_3e <- function(explor = TRUE) {
  w <- sz("gsea_nes")
  g <- fig_gsea_nes(layout = "dodged")
  if (!is.null(g)) save_plot(g, "Fig3E_gsea_nes_by_subset", w[1], w[2])
  # The three-facet alternative is kept as an exploratory output so the two can
  # be compared again if the layout is ever reopened. It is not the panel, so
  # explor = FALSE skips it when only the panel is being redrawn.
  if (explor) {
    g <- fig_gsea_nes(layout = "facet")
    if (!is.null(g)) save_plot(g, "explor_gsea_nes_facet", 7.0, 2.8)
  }
  invisible(TRUE)
}

panel_s1d <- function() {
  w <- named_gene_size()
  g <- fig_named_gene_violins()
  if (!is.null(g)) save_plot(g, "FigS1D_named_genes_by_subset", w[1], w[2])
  invisible(TRUE)
}

# The two exploratory Figure S1 outputs. Not manuscript panels, which is why
# they are one function rather than two and why their names start explor_.
panels_s1_exploratory <- function() {
  grp <- fib_grp()
  lv <- levels(mouse_fib@meta.data[[grp]])
  d <- expr_long(mouse_fib, VLN_GENES_MOUSE, grp) %>%
    mutate(group = factor(group, levels = lv))
  w <- sz("violin")
  save_plot(violin_panel(d, "group", "expression", "gene", ncol = 4,
                         ylab = "Normalised expression") +
              theme(strip.text = element_text(face = "italic")),
            "explor_violin_subset_markers_s1", w[1], w[2])

  w <- sz("fib_features_s1")
  save_plot(feature_plot(mouse_fib, FIG$fib_feature_genes_s1, group_by = grp,
                         ncol = 3, width_in = w[1]),
            "explor_featureplots_fibroblasts_s1", w[1], w[2])
  invisible(TRUE)
}

# --- by name ------------------------------------------------------------------
# The lookup is defined AFTER the functions, so every value here is a function
# object that already exists; a typo in a name is an error at load, not a
# confusing failure half way through a run.
# RENUMBERED 2026-10-08. The GSEA panel moved out of the supplement into
# Figure 3 and pushed the module scores down a letter:
#     Fig S1D (GSEA)          -> Fig 3E
#     Fig 3E  (module scores) -> Fig 3F
#     Fig S1E (named genes)   -> Fig S1D
# Panel letters belong in THESE names and in the output file names, and
# nowhere else. The constants that drive the named-gene panel were called
# S1E_* and had to be renamed with it, which is why they are now
# NAMED_GENE_* -- matching fig_named_gene_violins() and
# PANEL_SIZE$named_gene_vln, which never carried a letter and never needed
# touching. Keep it that way: a letter in a content name is a renumbering bill.
PANEL_INDEX <- list(
  "3A1" = panel_3a1, "3A2" = panel_3a2, "3B" = panel_3b, "3C" = panel_3c,
  "3D"  = panel_3d,  "3E"  = panel_3e,  "3F"  = panel_3f,
  "S1A" = panel_s1a, "S1B" = panel_s1b, "S1C" = panel_s1c,
  "S1D" = panel_s1d)

list_panels <- function() {
  message("panels(\"3B\") or panel_3b() -- one panel, one file:")
  for (n in names(PANEL_INDEX)) message("  ", n)
  message("fig_panels() draws all of them plus the exploratory outputs.")
  invisible(names(PANEL_INDEX))
}

# panels("3B", "S1D") -- redraw a few without the rest. Case does not matter
# and "Fig" is optional, so panels("figS1c") works.
panels <- function(...) {
  want <- toupper(gsub("^FIG", "", toupper(unlist(c(...)))))
  bad <- setdiff(want, names(PANEL_INDEX))
  if (length(bad))
    stop("no such panel: ", paste(bad, collapse = ", "),
         "\nknown: ", paste(names(PANEL_INDEX), collapse = ", "))
  for (n in want) { message(n); PANEL_INDEX[[n]]() }
  invisible(TRUE)
}

fig_panels <- function() {
  ensure("mouse_all"); ensure("mouse_fib")
  message("figures_final.R ", FIG_VERSION,
          " -- FINAL manuscript panels; Fig3 A/B/C/D/E = UMAPs / dot plot / ",
          "feature plots / GO / module scores")
  message("  from : ", tryCatch(normalizePath("figures.R"),
                                error = function(e) "?"))
  message("  to   : ", normalizePath(FIGDIR, mustWork = FALSE))
  message("  Fig3C: ", FIG$fib_feature_ncol, " columns, labels \"",
          FIG$feature_label, "\", panel ",
          paste(sz("fib_features"), collapse = " x "), " in")

  # --- Figure 3 ---------------------------------------------------------------
  # 3C first, because 3A2 and 3A1 are matched to its dot size. That is no
  # longer a requirement -- fib_dot_target() recomputes the geometry for
  # whoever asks -- but drawing them in this order keeps the log readable.
  panel_3c(); panel_3a2(); panel_3a1(); panel_3b(); panel_3d()
  panel_3e(); panel_3f()

  # --- Figure S1 --------------------------------------------------------------
  panel_s1a(); panel_s1b(); panel_s1c(); panel_s1d()
  panels_s1_exploratory()

  message("Done. Place every panel at 100 % -- these are final sizes.")
  invisible(TRUE)
}
