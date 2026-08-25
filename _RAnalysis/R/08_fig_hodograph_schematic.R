# 08_fig_hodograph_schematic.R — Fig 1: canonical CoM velocity hodographs. A methods
# schematic that introduces how to read a hodograph and its ROTATION SENSE before the
# data figures use it (Figs 3, 5 and S4). Two idealized loops in the (fore-aft velocity
# fluctuation, vertical velocity) plane: a pendular walk-type loop (counter-clockwise,
# positive enclosed area, kinetic and potential energy out of phase) and a bouncing
# run-type loop (clockwise, negative area, kinetic and potential energy in phase).
# Grounded running and aerial running are both bouncing (clockwise); see caption.
# Idealized, so the axes carry no numeric scale.
suppressPackageStartupMessages({library(dplyr); library(ggplot2)})
if (!exists("theme_daley")) source("R/lib_theme.R")

th <- seq(0, 1, length.out = 240)
# Start each loop at the CANONICAL touchdown phase for its gait rather than a shared point:
#   walk (pendular): touchdown at fastest forward speed with vertical velocity ~ 0, then the
#     CoM vaults UP -> dot on the right (3 o'clock), loop runs counter-clockwise.
#   run (bouncing): touchdown while the CoM is still descending from flight, dropping to
#     maximum downward velocity before rebounding -> dot at lower-right (~4:30), loop clockwise.
# The start angles match where touchdown falls in the real steady hodographs (Fig 5): walk
# vertical velocity ~ 0 at contact; aerial-run contact near -55 deg (lower-right, descending).
mk <- function(sense, lab, phi0) {
  s <- if (sense == "CCW") 1 else -1                 # +1 CCW (angle increases), -1 CW
  a <- phi0 + s * 2*pi*th
  data.frame(phase = th*100, panel = lab, x = 0.85 * cos(a), y = 1.1 * sin(a))
}
LAB_CCW <- "Pendular (walk-type)\nCCW, area > 0\nKE, PE out of phase"
LAB_CW  <- "Bouncing (run-type)\nCW, area < 0\nKE, PE in phase"
loops <- rbind(mk("CCW", LAB_CCW, 0), mk("CW", LAB_CW, -pi/4))
loops$panel <- factor(loops$panel, levels = c(LAB_CCW, LAB_CW))
pal   <- c("#3D4A5C","#1C9DA8","#C0398B","#D9A23B")
start <- loops %>% group_by(panel) %>% slice_min(phase, n = 1) %>% ungroup()
arrw  <- loops %>% group_by(panel) %>% filter(phase >= 6, phase <= 14) %>% ungroup()  # direction, just after touchdown

p <- ggplot(loops, aes(x, y, colour = phase)) +
  geom_hline(yintercept = 0, linetype = 3, colour = "grey70") +
  geom_vline(xintercept = 0, linetype = 3, colour = "grey70") +
  geom_path(aes(group = panel), linewidth = 0.7) +
  geom_path(data = arrw, aes(group = panel), linewidth = 0.7,
            arrow = arrow(length = unit(1.6, "mm"), type = "closed")) +
  geom_point(data = start, colour = "black", size = 1.3) +
  scale_colour_gradientn(colours = pal, name = "cycle %", limits = c(0,100)) +
  facet_wrap(~panel, ncol = 2) + coord_equal() +
  labs(x = LAB_FA_VEL, y = LAB_VERT_VEL) +
  theme_bio() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        strip.text = element_text(size = BIO_BASE, face = "bold", lineheight = 1.0),
        panel.background = element_rect(fill = "transparent", colour = NA),
        legend.position = "bottom", legend.title = element_text(vjust = 1),
        legend.key.height = unit(6, "pt"), legend.key.width = unit(26, "pt"),
        legend.box.spacing = unit(2, "pt"), legend.margin = margin(0, 0, 0, 0),
        axis.title.x = element_text(margin = margin(t = 2)),
        # A strip label is clipped to its own panel, and these run the full panel width.
        strip.clip = "off", panel.spacing.x = unit(2, "pt")) +
  guides(colour = guide_colourbar(title.position = "left"))
# Sized to the content, not to a round number: coord_equal() holds the panel aspect, so surplus
# canvas in either direction becomes blank margin instead of larger panels. The height below is
# measured against the two panels at single-column width. See the same note in
# 14_fig_hodographs.R.
bio_save("Fig1_HodographSchematic", p, width = BIO_W1, height = 2.55)
cat("08_fig_hodograph_schematic: Fig1_HodographSchematic (canonical hodograph schematic).\n")
