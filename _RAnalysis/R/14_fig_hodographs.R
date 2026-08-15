# 14_fig_hodographs.R — Fig 5: steady-gait CoM velocity hodographs, averaged over
# the STRIDE cycle with Hilbert continuous-phase registration. Fore-aft velocity
# fluctuation (about each stride's mean) on x, vertical velocity on y, both
# /sqrt(gL0); loops coloured by percent of the STRIDE, 0 to 100. The Hilbert registration
# indexes a stride as two consecutive 0-100 step phases, so the phase variable runs 0 to 200 and
# is halved for display. Black dot = stride touchdown.
#
# Panels are gait x speed half, the SAME grouping Fig 7 uses (17_fig_mean_traces.R): each gait is
# split slow/fast at its median Froude over all classified observations, the threshold
# 05_mean_traces.R writes to output/speedbin_thresholds_R.csv. Reading the same split on the loops
# and on the force and energy traces shows how loop shape shifts with speed within a gait as well
# as between gaits. The six panels sit in one row ordered by median speed.
#
# Stride-cycle registration is used (rather than a single step) because the CoM returns to its
# state over a full stride, so the mean loop closes; plain percent-of-cycle step averaging over
# the heterogeneous multi-study sample leaves the loop open (phase jitter), and Hilbert
# registration aligns each stride by its vertical-velocity oscillation.
suppressPackageStartupMessages({library(dplyr); library(ggplot2)})
if (!exists("theme_daley")) source("R/lib_theme.R")
if (!exists("pr_registered_transitions_hilbert")) source("R/lib_phase_register.R")
ctf <- "data/cycleTracesStep.csv"
thf <- "output/speedbin_thresholds_R.csv"
if (!file.exists(ctf)) { cat("14_fig_hodographs: cycleTracesStep.csv not found.\n") } else {
  steadyKeys <- c("walk","groundedRun","aerialRun")
  steadyLab  <- c(walk="Walk", groundedRun="Grounded run", aerialRun="Aerial run")
  SORD <- c("slow","fast")
  ctStep <- readr::read_csv(ctf, show_col_types = FALSE)

  # the gait x speed split Fig 7 uses, read from the file 05_mean_traces.R writes so the two
  # figures cannot drift apart
  spd_thresh <- if (file.exists(thf)) readr::read_csv(thf, show_col_types = FALSE) else
    step %>% filter(!is.na(gaitObjective)) %>% group_by(gait = as.character(gaitObjective)) %>%
      summarise(medianFroude = median(Froude, na.rm = TRUE), .groups = "drop")

  # steady strides labelled by gait and speed half, registered over the stride by Hilbert phase
  # (0-200, halved to percent of stride for the colour scale)
  steadyStr <- stride %>%
    filter(as.character(gaitObjective) %in% steadyKeys, accClass == "steady") %>%
    mutate(gait = as.character(gaitObjective)) %>%
    left_join(spd_thresh, by = "gait") %>%
    mutate(spd = ifelse(Froude >= medianFroude, "fast", "slow")) %>%
    transmute(boutID, strideIndex, transClass = paste(gait, spd, sep = "|"))
  reg <- pr_registered_transitions_hilbert(ctStep, steadyStr)
  nstr <- reg %>% distinct(strideID, transClass) %>% count(transClass)
  meanW <- reg %>% group_by(transClass, phase) %>%
    summarise(vfa_c = mean(vfa_c, na.rm = TRUE), vvert_n = mean(vvert_n, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(gait = sub("\\|.*$", "", transClass), spd = sub("^.*\\|", "", transClass),
           pctStride = phase / 2)   # phase spans two step cycles = one stride

  # observed rotation sense and stride count, per gait (quoted in Results) and per gait x
  # speed cell (the panels of the figure)
  steadyRot <- stride %>%
    filter(as.character(gaitObjective) %in% steadyKeys, accClass == "steady",
           is.finite(hodoArea_n)) %>%
    mutate(gait = as.character(gaitObjective)) %>%
    left_join(spd_thresh, by = "gait") %>%
    mutate(spd = ifelse(Froude >= medianFroude, "fast", "slow"))
  # The median SIGNED area is emitted per cell, not just the rotation sense. Ordered by speed it
  # falls monotonically through zero, and the walking loop therefore SHRINKS with speed while the
  # running loops grow: the area is heading for the crossover, not away from it.
  summarise_rot <- function(d) d %>%
    summarise(n_strides = dplyr::n(), pct_ccw = round(100 * mean(hodoArea_n > 0)),
              median_u = round(median(meanSpeed_n), 2),
              median_signed_area = round(median(hodoArea_n), 4), .groups = "drop") %>%
    mutate(sense = ifelse(pct_ccw >= 50, "CCW", "CW"),
           pct_majority = ifelse(pct_ccw >= 50, pct_ccw, 100 - pct_ccw))
  rotObs <- summarise_rot(steadyRot %>% group_by(gait, spd))
  rotGait <- summarise_rot(steadyRot %>% group_by(gait)) %>% mutate(spd = "all", .after = gait)
  readr::write_csv(bind_rows(rotGait, rotObs), "output/fig4_rotation_sense_observed.csv")

  # Panel label: gait, speed half and the observed majority rotation sense with the percentage
  # of strides turning that way, all read off the figure so the legend carries none of it.
  #
  # One row of six, ordered by the cell's median dimensionless speed, so the panels run from the
  # slowest walk to the fastest aerial run and the change in loop shape with speed reads left to
  # right across gaits as well as within them.
  key_of  <- paste(rotObs$gait, rotObs$spd, sep = "|")
  cellSns <- setNames(rotObs$sense, key_of)
  cellPct <- setNames(rotObs$pct_majority, key_of)
  celllab <- function(key) {
    g <- sub("\\|.*$", "", key); sp <- sub("^.*\\|", "", key)
    sprintf("%s, %s\n%d%% %s",
            steadyLab[g], sp, cellPct[key], cellSns[key])
  }
  cellOrder <- key_of[order(rotObs$median_u)]
  meanW$facet <- factor(celllab(meanW$transClass), levels = celllab(cellOrder))
  start <- meanW %>% group_by(facet) %>% slice_min(phase, n = 1) %>% ungroup()

  p4 <- ggplot(meanW, aes(vfa_c, vvert_n, colour = pctStride, group = facet)) +
    geom_hline(yintercept = 0, linetype = 3, colour = "grey70") +
    geom_vline(xintercept = 0, linetype = 3, colour = "grey70") +
    geom_path(linewidth = 1.2) +
    geom_point(data = start, aes(vfa_c, vvert_n), colour = "black", size = 2, inherit.aes = FALSE) +
    scale_colour_gradientn(colours = c("#3D4A5C","#1C9DA8","#C0398B","#D9A23B"),
                           name = "stride %", limits = c(0, 100)) +
    scale_x_continuous(breaks = breaks3) + scale_y_continuous(breaks = breaks3) +
    facet_wrap(~facet, nrow = 1) + coord_equal() +
    labs(x = LAB_FA_VEL, y = LAB_VERT_VEL) +
    theme_daley() +
    theme(strip.text = element_text(size = 8, lineheight = 1.0),
          panel.background = element_rect(fill = "transparent", colour = NA),
          strip.background = element_rect(fill = "transparent", colour = NA),
          panel.spacing = grid::unit(8, "pt"))
  # coord_equal() fixes the panel aspect, so the canvas has to be sized to what the panels
  # actually occupy: a wider canvas pads them with blank rather than enlarging them, and the page
  # then scales the whole figure down. Re-measure the content bounding box if the panel count, the
  # legend or the axis labels change.
  ggsave("output/Fig5_SteadyHodographs.pdf", p4, width = 8.7, height = 4.7, device = cairo_pdf)
  ggsave("output/Fig5_SteadyHodographs.png", p4, width = 8.7, height = 4.7, dpi = 300, bg = "white")
  cat(sprintf("14_fig_hodographs: Fig5_SteadyHodographs, stride-Hilbert, gait x speed (%s).\n",
              paste(sprintf("%s=%d", nstr$transClass, nstr$n), collapse = ", ")))
  print(as.data.frame(rotObs))
}
