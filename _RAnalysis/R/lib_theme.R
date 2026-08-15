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

theme_daley <- function(base_size = 11) {
  theme_classic(base_size = base_size, base_family = .base_family) +
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
LAB_STEPLEN  <- expression("step length / " * italic(L)[0])
LAB_STEPFREQ <- expression("step frequency  " * sqrt(italic(L)[0] / italic(g)))

## 3 axis ticks including zero (for the crowded hodograph axes): symmetric about 0.
breaks3 <- function(lims) { m <- max(abs(lims), na.rm = TRUE); b <- signif(0.6 * m, 1); c(-b, 0, b) }
