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
LAB_CCW <- "Pendular (walk-type)\nCCW, loop area > 0, KE & PE out of phase"
LAB_CW  <- "Bouncing (run-type)\nCW, loop area < 0, KE & PE in phase"
loops <- rbind(mk("CCW", LAB_CCW, 0), mk("CW", LAB_CW, -pi/4))
loops$panel <- factor(loops$panel, levels = c(LAB_CCW, LAB_CW))
pal   <- c("#3D4A5C","#1C9DA8","#C0398B","#D9A23B")
start <- loops %>% group_by(panel) %>% slice_min(phase, n = 1) %>% ungroup()
arrw  <- loops %>% group_by(panel) %>% filter(phase >= 6, phase <= 14) %>% ungroup()  # direction, just after touchdown

p <- ggplot(loops, aes(x, y, colour = phase)) +
  geom_hline(yintercept = 0, linetype = 3, colour = "grey70") +
  geom_vline(xintercept = 0, linetype = 3, colour = "grey70") +
  geom_path(aes(group = panel), linewidth = 1.2) +
  geom_path(data = arrw, aes(group = panel), linewidth = 1.2,
            arrow = arrow(length = unit(3, "mm"), type = "closed")) +
  geom_point(data = start, colour = "black", size = 2.4) +
  scale_colour_gradientn(colours = pal, name = "cycle %", limits = c(0,100)) +
  facet_wrap(~panel, ncol = 2) + coord_equal() +
  labs(x = LAB_FA_VEL, y = LAB_VERT_VEL) +
  theme_daley() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        strip.text = element_text(size = 9, lineheight = 0.95),
        panel.background = element_rect(fill = "transparent", colour = NA))
# Sized to the content, not to a round number: coord_equal() holds the panel aspect, so extra
# canvas width becomes blank margin rather than larger panels. See the same note in
# 14_fig_hodographs.R.
ggsave("output/Fig1_HodographSchematic.pdf", p, width = 7.0, height = 4.6, device = cairo_pdf)
ggsave("output/Fig1_HodographSchematic.png", p, width = 7.0, height = 4.6, dpi = 300, bg = "white")
cat("08_fig_hodograph_schematic: Fig1_HodographSchematic (canonical hodograph schematic).\n")
