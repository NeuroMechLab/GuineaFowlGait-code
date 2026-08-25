# lib_theme.R — shared plotting layer (Neuromech Lab house style)
# Colourblind-safe palette (Wong/Tol), no gridlines / box / mirrored axes,
# separate panels (no overlaid twin axes). Redundant colour + shape + linetype.

suppressPackageStartupMessages({library(ggplot2)})

# Gait palette: the three gaits, classified by two features (aerial phase; hodograph
# rotation sense). walk (navy, pendular/vaulting), grounded run (teal, bouncing, no
# flight), aerial run (magenta, bouncing with flight). Rotation-sense
# (pendular/bouncing) palette provided for the gait-space figure.
GAIT_LEVELS <- c("walk","groundedRun","aerialRun")
GAIT_LABELS <- c(walk="Walk", groundedRun="Grounded run", aerialRun="Aerial run")
GAIT_COLORS <- c(walk="#3D4A5C", groundedRun="#1C9DA8", aerialRun="#C0398B")
GAIT_SHAPES <- c(walk=16, groundedRun=17, aerialRun=15)
# Redundant encoding for outlines is LINE THICKNESS, never a dash pattern: a dashed or dotted
# outline breaks up a histogram edge and reads as a different kind of data. Thickest for walking,
# which is the smallest sample and the one most easily hidden under the others.
GAIT_LWD    <- c(walk=1.15, groundedRun=0.75, aerialRun=0.45)
ROT_COLORS  <- c(`pendular (walk)`="#3D4A5C", `bouncing (run)`="#C0398B")

# Four-cell gait class (02_clean.R): the same two features crossed, with the pendular
# flight steps split out of aerialRun as their own cell so it can be quantified rather
# than asserted empty. Used by the four-cell summary tables (20_table_gait4_summary.R).
# The pendular run takes the orange of the Wong palette, distinct from the three gait colours.
# The cell is named for the gait Srinivasan & Ruina (2006) predicted, which is the name the
# paper, the figures, the tables and the archived data all use.
GAIT4_LEVELS <- c("walkGrounded","groundedRun","aerialRun","pendularRun")
GAIT4_LABELS <- c(walkGrounded="Walk (grounded)", groundedRun="Grounded run",
                  aerialRun="Aerial run", pendularRun="Pendular run")
GAIT4_COLORS <- c(walkGrounded="#3D4A5C", groundedRun="#1C9DA8",
                  aerialRun="#C0398B", pendularRun="#E69F00")
GAIT4_SHAPES <- c(walkGrounded=16, groundedRun=17, aerialRun=15, pendularRun=4)

.have_roboto <- "Roboto" %in% tryCatch(systemfonts::system_fonts()$family, error=function(e) character(0))
.base_family <- if (isTRUE(.have_roboto)) "Roboto" else ""

theme_daley <- function(base_size = 11, base_family = .base_family) {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line   = element_line(colour = "grey20", linewidth = 0.4),
      axis.ticks  = element_line(colour = "grey20", linewidth = 0.4),
      panel.grid  = element_blank(),
      panel.border= element_blank(),
      strip.background = element_blank(),
      strip.text  = element_text(face = "bold"),
      legend.key  = element_blank(),
      legend.title= element_text(face = "bold"),
      plot.title  = element_text(face = "bold")
    )
}

scale_color_gait <- function(...) scale_color_manual(values = GAIT_COLORS, labels = GAIT_LABELS,
                                                     breaks = GAIT_LEVELS, name = "Gait", ...)
scale_fill_gait  <- function(...) scale_fill_manual(values = GAIT_COLORS, labels = GAIT_LABELS,
                                                    breaks = GAIT_LEVELS, name = "Gait", ...)
scale_shape_gait <- function(...) scale_shape_manual(values = GAIT_SHAPES, labels = GAIT_LABELS,
                                                     breaks = GAIT_LEVELS, name = "Gait", ...)

## Normalization axis labels — one consistent format everywhere (proper square-root,
## italic g, italic L with subscript 0, italic m). Use these plotmath objects in
## labs()/xlab()/ylab() instead of plain text like "(/sqrt(gL0))" or "/mgL0".
LAB_U        <- expression("dimensionless speed  " * italic(v) / sqrt(italic(g) * italic(L)[0]))
LAB_FA_VEL   <- expression("fore-aft velocity fluctuation / " * sqrt(italic(g) * italic(L)[0]))
LAB_VERT_VEL <- expression("vertical velocity / " * sqrt(italic(g) * italic(L)[0]))
LAB_HODOAREA <- expression("hodograph signed area / " * italic(g) * italic(L)[0])
LAB_ENERGY   <- expression("energy fluctuation / " * italic(m) * italic(g) * italic(L)[0])
## Two-line form, for a panel shorter than the one-line label is long. A y title taller than its
## panel pushes into the panel above it once stacked rows are drawn close together.
LAB_ENERGY_2L <- expression(atop("energy fluctuation", "/ " * italic(m) * italic(g) * italic(L)[0]))
LAB_STEPLEN  <- expression("step length / " * italic(L)[0])
LAB_STEPFREQ <- expression("step frequency  " * sqrt(italic(L)[0] / italic(g)))

## 3 axis ticks including zero (for the crowded hodograph axes): symmetric about 0.
breaks3 <- function(lims) { m <- max(abs(lims), na.rm = TRUE); b <- signif(0.6 * m, 1); c(-b, 0, b) }
## The same pair without the zero tick, for an axis narrower than three labels: the dotted zero
## line already marks the origin.
breaks2 <- function(lims) { m <- max(abs(lims), na.rm = TRUE); b <- signif(0.6 * m, 1); c(-b, b) }


## ---- Biology Open figure specification --------------------------------------------------
## Widths and the height bound are the journal's; the type sizes are its 8 pt body labelling
## and 12 pt bold capital panel tags. Export each figure at the width it is printed at, so a
## nominal point size is the printed point size: a canvas wider than the placed width shrinks
## every label by the ratio between the two.
BIO_W1   <- 88  / 25.4   # single-column width, inches
BIO_W2   <- 183 / 25.4   # double-column width, inches
BIO_HMAX <- 210 / 25.4   # maximum height, inches
BIO_BASE <- 8            # pt, every label other than the panel tags
BIO_TAG  <- 12           # pt, panel tags, bold capitals
BIO_FONT <- "Arial"

theme_bio <- function(base_size = BIO_BASE) {
  theme_daley(base_size = base_size, base_family = BIO_FONT) +
    theme(
      axis.text    = element_text(size = base_size, colour = "grey20"),
      axis.title   = element_text(size = base_size),
      legend.text  = element_text(size = base_size),
      legend.title = element_text(size = base_size, face = "bold"),
      strip.text   = element_text(size = base_size, face = "bold"),
      plot.tag     = element_text(family = BIO_FONT, face = "bold", size = BIO_TAG),
      plot.margin  = margin(2, 2, 2, 2)
    )
}

## Panels stacked on one shared x axis: the bottom panel carries the axis and the panels above
## drop their tick labels, title and ticks, which is what lets the rows sit close together.
## Apply only where the panels above and below are drawn on the same x range.
bio_drop_x <- function(bottom = 0) {
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.ticks.x = element_blank(), plot.margin = margin(2, 2, bottom, 2))
}

## Legend inside the panel, positioned in normalized panel coordinates.
bio_inset_legend <- function(x = 0.99, y = 0.02, just = c(1, 0), base_size = BIO_BASE) {
  theme(legend.position = "inside", legend.position.inside = c(x, y),
        legend.justification = just,
        legend.background = element_blank(), legend.key = element_blank(),
        legend.margin = margin(0, 0, 0, 0), legend.spacing.y = unit(1, "pt"),
        legend.key.size = unit(base_size + 2, "pt"),
        legend.text = element_text(size = base_size),
        legend.title = element_text(size = base_size, face = "bold"))
}

## Every figure is written as a vector PDF at the printed size and as a PNG for the built
## manuscript. cairo_pdf writes RGB and embeds no colour profile, which is what the journal asks
## for, and it embeds the Arial faces so the labels survive off this machine.
bio_save <- function(name, plot, width, height, dir = "output") {
  stopifnot(width  <= BIO_W2 + 1e-6, height <= BIO_HMAX + 1e-6)
  ggsave(file.path(dir, paste0(name, ".pdf")), plot, width = width, height = height,
         device = cairo_pdf, bg = "white")
  ggsave(file.path(dir, paste0(name, ".png")), plot, width = width, height = height,
         dpi = 300, bg = "white")
  cat(sprintf("  %s: %.1f x %.1f mm\n", name, width * 25.4, height * 25.4))
}
