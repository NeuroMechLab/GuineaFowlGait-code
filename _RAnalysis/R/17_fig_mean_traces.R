# 17_fig_mean_traces.R — Fig 7 (steady) plus the accel/decel supplements: mean
# STEP-cycle traces by GAIT and speed, one figure per steadiness class. Columns are
# gait x speed bin (walk / grounded run / aerial run, each split slow/fast at the
# gait's median Froude). Splitting within gait keeps each column gait-pure, so the
# flight phase survives the averaging, and gives a speed progression across columns.
# Traces are group means on a real-time axis with a Student-t CI, computed in R by
# 05_mean_traces.R (grouping classGaitSpd) so the averaging groups follow the gait label this
# pipeline assigns. Rows: vertical GRF, fore-aft GRF, CoM energy (KE, PE); zero lines on
# fore-aft & energy.
#
# The speed and sample size of each column sit in a FOURTH row beneath the energy panels, as a
# table of Froude, dimensionless speed, absolute speed (each mean [min max]) and n steps, one
# column per group. A full row gives the four lines of numbers the page width they need to stay
# legible at the printed size.
#
# The canvas is exported at the width the figure is printed at, so a nominal point size is the
# printed point size and the journal's 8 pt labelling is drawn at 8 pt. A canvas wider than the
# placed width would shrink every label by the ratio between the two.
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2); library(patchwork)})
if (!exists("theme_daley")) source("R/lib_theme.R")
if (!exists("meanTraces")) source("R/05_mean_traces.R")
if (is.null(meanTraces)) { cat("17_fig_mean_traces: no mean traces available.\n") } else {
  MT <- meanTraces %>%
    filter(grouping == "classGaitSpd", cycleType == "step") %>%
    separate(frbin, into = c("gait","spd"), sep = "\\|", remove = FALSE)
  glab <- c(walk="Walk", groundedRun="Grounded run", aerialRun="Aerial run")
  GORD <- names(glab); SORD <- c("slow","fast")
  # column order: gait (walk<grounded<aerial) then speed (slow<fast)
  colLevels <- as.vector(t(outer(GORD, SORD, paste, sep="|")))

  # Page geometry (see the header note): double-column width, and the journal's body type size.
  FIG_W <- BIO_W2; FIG_H <- 5.0     # export canvas, inches
  BASE  <- BIO_BASE

  # 05_mean_traces.R writes the per-cell summary over the EXACT steps behind each mean trace,
  # so the column label, the n and the Fr/u/v ranges describe the drawn curve rather than an
  # approximation. The strip carries the gait and speed group; the numbers go to the table row.
  # Each column carries its own free x scale, so the last tick of one panel and the first tick
  # of the next sit side by side: at this type size they collide ("200" and "0" reading as
  # "2000") unless the panels are held apart.
  GAP  <- theme(panel.spacing.x = grid::unit(0.9, "lines"))
  # Four ticks, not the default six. Each column spans a different cycle duration, so ggplot
  # picks the break count per panel from its own range; six labels do not fit the width of one
  # panel of six at this type size, and run together ("0 50100150200250").
  XBRK <- scale_x_continuous(n.breaks = 4)
  TAB_ROWS <- c("Fr", "u", "v (m/s)", "n steps")
  ann <- meanTracesAnn %>% filter(grouping == "classGaitSpd") %>%
    mutate(col = paste(gait, spd, sep = "|"),
           # Two lines, not one: "Grounded run (slow)" set at BASE pt is wider than one panel
           # of a six-column grid, and ggplot clips a strip label rather than shrinking it.
           lab = sprintf("%s\n(%s)", glab[as.character(gait)], spd)) %>%
    transmute(accClass = grp, col, lab,
              `Fr`        = sprintf("%.2f [%.2f %.2f]", Fr_mean, Fr_min, Fr_max),
              `u`         = sprintf("%.2f [%.2f %.2f]", u_mean, u_min, u_max),
              `v (m/s)`   = sprintf("%.2f [%.2f %.2f]", v_mean, v_min, v_max),
              `n steps`   = sprintf("%d", n))
  anntab <- ann %>%
    tidyr::pivot_longer(all_of(TAB_ROWS), names_to = "row", values_to = "val") %>%
    mutate(ypos = length(TAB_ROWS) - match(row, TAB_ROWS) + 1)

  # The KE and PE key sits in whichever corner of the energy row that figure leaves empty:
  # bottom right for the steady and accelerating traces, top right for the decelerating ones,
  # whose kinetic-energy trace falls through the bottom right of the fastest column.
  make_fig <- function(cls, outname, key_x = 0.99, key_y = 0.02, key_just = c(1, 0)) {
    d <- MT %>% filter(grp == cls)
    if (!nrow(d)) { cat(sprintf("  (no traces for %s)\n", cls)); return(invisible()) }
    cols <- intersect(colLevels, unique(d$frbin))
    labmap <- setNames(ann$lab[match(paste(cls,cols), paste(ann$accClass, ann$col))], cols)
    labmap[is.na(labmap)] <- cols[is.na(labmap)]
    d$frbin <- factor(d$frbin, levels = cols)
    panel <- function(ch, ylab, col, zero=FALSE, strip=FALSE) {
      dd <- d %>% filter(channel == ch)
      p <- ggplot(dd, aes(t_ms)) +
        geom_ribbon(aes(ymin = lo, ymax = hi), fill = col, alpha = 0.25) +
        geom_line(aes(y = mean), colour = col, linewidth = 0.7)
      if (zero) p <- p + geom_hline(yintercept = 0, linetype = 3, colour = "grey55")
      p + XBRK + facet_grid(. ~ frbin, scales = "free_x",
                     labeller = if (strip) labeller(frbin = labmap) else label_value) +
        labs(x = NULL, y = ylab) + theme_bio(BASE) +
        theme(strip.text = if (strip) element_text(size = BASE, face = "bold") else element_blank(),
              strip.background = element_blank()) +
        GAP
    }
    pFz  <- panel("Fz_BW",  "vertical GRF (BW)", "#3D4A5C", strip=TRUE)
    pFfa <- panel("Ffa_BW", "fore-aft GRF (BW)", "#0B6E4F", zero=TRUE)
    en <- d %>% filter(channel %in% c("KE_n","PE_n")) %>% mutate(E = sub("_n","",channel))
    pEn <- ggplot(en, aes(t_ms, colour = E, fill = E)) +
      geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, colour = NA) +
      geom_line(aes(y = mean), linewidth = 0.7) + geom_hline(yintercept = 0, linetype = 3, colour = "grey55") +
      scale_colour_manual(values = c(KE="#0077BB", PE="#C0398B"), name=NULL) +
      scale_fill_manual(values = c(KE="#0077BB", PE="#C0398B"), name=NULL) +
      XBRK + facet_grid(. ~ frbin, scales = "free_x") +
      labs(x = "time (ms)", y = LAB_ENERGY_2L) + theme_bio(BASE) +
      theme(strip.text = element_blank(), strip.background = element_blank()) + GAP +
      bio_inset_legend(key_x, key_y, just = key_just, base_size = BASE)

    # The speed-and-sample-size table. It is a facet_grid over the same column factor as the
    # trace panels, so patchwork aligns each cell under the group it describes; the row names
    # ride on the y axis, which puts them once at the left rather than in every cell.
    td <- anntab %>% filter(accClass == cls, col %in% cols) %>%
      mutate(frbin = factor(col, levels = cols))
    pTab <- ggplot(td, aes(x = 0, y = ypos, label = val)) +
      geom_hline(yintercept = length(TAB_ROWS) + 0.55, colour = "grey80", linewidth = 0.4) +
      geom_text(size = BASE * 0.352777778, family = BIO_FONT) +
      facet_grid(. ~ frbin) +
      scale_x_continuous(limits = c(-1, 1), expand = c(0, 0)) +
      scale_y_continuous(breaks = seq_along(TAB_ROWS), labels = rev(TAB_ROWS),
                         limits = c(0.4, length(TAB_ROWS) + 0.7), expand = c(0, 0)) +
      labs(x = NULL, y = NULL) + theme_bio(BASE) +
      theme(strip.text = element_blank(), axis.line = element_blank(),
            axis.ticks = element_blank(), axis.text.x = element_blank()) + GAP

    # Every column draws its three rows on one time range, so the axis is carried by the energy
    # row and the rows above it sit together. The KE and PE key sits inside that row.
    fig <- (pFz + bio_drop_x()) / (pFfa + bio_drop_x()) / pEn / pTab +
      plot_layout(heights = c(1, 1, 1, 0.55))
    bio_save(outname, fig, width = FIG_W, height = FIG_H)
    cat(sprintf("    columns %s\n", paste(cols, collapse=", ")))
  }
  cat("11_fig_meantraces (by gait x speed bin):\n")
  make_fig("steady",       "Fig7_MeanTraces_Steady")
  make_fig("accelerating", "SFig_MeanTraces_Accel")
  make_fig("decelerating", "SFig_MeanTraces_Decel", key_y = 0.98, key_just = c(1, 1))
}
