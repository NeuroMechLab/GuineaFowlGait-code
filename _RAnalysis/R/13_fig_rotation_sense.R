# 13_fig_rotation_sense.R — Fig 3: hodograph ROTATION DIRECTION as an alternative
# gait discriminator. Over a step the CoM velocity vector traces a loop in the
# (fore-aft fluctuation, vertical) plane; the SIGN of the loop's enclosed area is
# its sense of rotation. Walking-type (pendular vaulting) and running-type
# (bouncing) loops turn opposite ways, so the sign is a binary walk-vs-run label
# that is independent of the energy-phase (recovery/congruity) classification.
# We compute the signed area per step, tabulate the rotation sense by gait, and fit a
# flexible model for the reversal: a binomial GAM P(running-type) ~ s(dimensionless speed) +
# s(individual, random effect). REML estimates the smoothing penalty, so the data set the
# steepness of the walk-to-run reversal; the fitted effective degrees of freedom are recorded
# with the model note.
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2); library(patchwork)})
if (!exists("theme_daley")) source("R/lib_theme.R")
if (!exists("signed_area")) source("R/lib_hodograph.R")
ctf <- "data/cycleTracesStep.csv"
if (!file.exists(ctf)) { cat("13_fig_rotation_sense: cycleTracesStep.csv not found.\n") } else {
  # The reported signed area is the production column hodoArea_n, computed once in MATLAB and
  # used by Table S1, the transition analysis and the feature matrices, so every part of the
  # paper reads the same number. The area is ALSO recomputed here from the exported cycle
  # traces, because the within-loop reversal diagnostics below need the per-sample trace
  # anyway; that recomputation serves as an independent check on the production value and its
  # agreement is emitted rather than assumed.
  sa <- readr::read_csv(ctf, show_col_types = FALSE) %>%
    mutate(stepIndex = strideIndex + stepInStride - 1) %>%
    inner_join(step %>% select(boutID, stepIndex, bird, subjectID, gaitObjective, accClass,
                               meanSpeed_n, Froude, hodoArea_n),
               by = c("boutID","stepIndex")) %>%
    filter(gaitObjective %in% GAIT_LEVELS) %>%
    group_by(boutID, stepIndex, bird, subjectID, gaitObjective, accClass, meanSpeed_n, Froude,
             hodoArea_n) %>%
    summarise(area_recomp = signed_area(vfa_n - mean(vfa_n), vvert_n), .groups = "drop") %>%
    mutate(gaitObjective = factor(gaitObjective, levels = GAIT_LEVELS),
           area = hodoArea_n,
           dir = ifelse(area > 0, "CCW (walking-type)", "CW (running-type)"),
           runType = as.integer(area < 0))
  readr::write_csv(data.frame(
    quantity = c("steps compared",
                 "max |cycle-trace recomputed area - production hodoArea_n|",
                 "Pearson r, recomputed vs production",
                 "steps where the two disagree in sign"),
    value = c(nrow(sa),
              signif(max(abs(sa$area_recomp - sa$hodoArea_n)), 3),
              round(cor(sa$area_recomp, sa$hodoArea_n), 6),
              sum(sign(sa$area_recomp) != sign(sa$hodoArea_n)))),
    "output_internal/hodoArea_recomputation_check.csv")

  # ---- how the signed area handles within-loop reversals and crossovers ----
  # The shoelace signed area is the NET circulation of the closed curve and is
  # translation invariant, so for a self-intersecting (figure-of-eight) loop it returns the
  # algebraic sum of the lobe areas: opposite-sense lobes cancel and the sign reports which
  # sense dominates by area. A local reversal therefore does not create a third category; it
  # only matters if it dominates. Three measures of the same curve, all from lib_hodograph.R:
  #   retro_frac    the fraction of the radius vector's turning angle running against the net
  #                 sense, which grades a reversal instead of categorising it
  #   n_crossings   how many times the TRAVERSED velocity path crosses itself, so a crossover is
  #                 counted as the property of the drawn curve that it is, not inferred from the
  #                 turning angle. The closing chord of the shoelace polygon is excluded, because
  #                 a step does not return to its starting velocity and the chord is a straight
  #                 line the bird never traversed (lib_hodograph.R); n_crossings_closed keeps the
  #                 polygon count so the size of that difference stays on record
  #   dom_share     the share of the loop's gross area lying in its net sense, which is the
  #                 signed area's own segment contributions split by sign, so a crossover can
  #                 be reported alongside how decisively the loop still turns one way
  rev <- readr::read_csv(ctf, show_col_types = FALSE) %>%
    mutate(stepIndex = strideIndex + stepInStride - 1) %>%
    inner_join(sa %>% select(boutID, stepIndex, gaitObjective, area), by = c("boutID","stepIndex")) %>%
    group_by(boutID, stepIndex, gaitObjective, area) %>%
    arrange(pct, .by_group = TRUE) %>%
    summarise(retro_frac  = loop_retro_frac(vfa_n, vvert_n),
              n_crossings = loop_self_crossings(vfa_n, vvert_n),
              n_crossings_closed = loop_self_crossings(vfa_n, vvert_n, closed = TRUE),
              dom_share   = dominant_area_share(vfa_n, vvert_n),
              closure_gap = sqrt((dplyr::last(vfa_n) - dplyr::first(vfa_n))^2 +
                                 (dplyr::last(vvert_n) - dplyr::first(vvert_n))^2) /
                            (diff(range(vfa_n)) + diff(range(vvert_n))),
              .groups = "drop") %>%
    mutate(crosses = n_crossings > 0)
  rev_tab <- rev %>% group_by(gait = gaitObjective) %>%
    summarise(n = n(),
              pct_with_reversal = round(100*mean(retro_frac > 0.01, na.rm = TRUE)),
              median_retro_pct  = round(100*median(retro_frac, na.rm = TRUE), 1),
              pct_retro_over_25 = round(100*mean(retro_frac > 0.25, na.rm = TRUE)),
              pct_with_crossover = round(100*mean(crosses, na.rm = TRUE)),
              median_dom_share_crossing = round(median(dom_share[crosses], na.rm = TRUE), 1),
              .groups = "drop")
  readr::write_csv(rev_tab, "output/hodograph_reversal_diagnostics.csv")
  cat("13_fig_rotation_sense: within-loop reversal and crossover diagnostics:\n")
  print(as.data.frame(rev_tab))

  # Overall figures quoted in Results and in the Discussion, emitted so every number in those
  # sentences is reproducible from the pipeline.
  rev_overall <- data.frame(
    n_steps                    = nrow(rev),
    pct_with_reversal          = round(100*mean(rev$retro_frac > 0.01, na.rm = TRUE)),
    median_retro_pct           = round(100*median(rev$retro_frac, na.rm = TRUE), 1),
    n_with_crossover           = sum(rev$crosses, na.rm = TRUE),
    pct_with_crossover         = round(100*mean(rev$crosses, na.rm = TRUE)),
    median_crossings_crossing  = stats::median(rev$n_crossings[rev$crosses], na.rm = TRUE),
    median_dom_share_crossing  = round(median(rev$dom_share[rev$crosses], na.rm = TRUE), 1),
    median_dom_share_simple    = round(median(rev$dom_share[!rev$crosses], na.rm = TRUE), 1),
    # the closing-chord contrast, so the choice of the traversed path is auditable
    n_with_crossover_closed    = sum(rev$n_crossings_closed > 0, na.rm = TRUE),
    median_closure_gap         = round(median(rev$closure_gap, na.rm = TRUE), 3))
  readr::write_csv(rev_overall, "output/hodograph_reversal_overall.csv")
  cat("  overall: "); print(as.data.frame(rev_overall))

  # Validation of the crossing test, for human review rather than for the paper: draw the loops
  # the test flags and the loops it does not, so the count can be checked by eye. The steps are
  # taken at evenly spaced ranks of the dominant-area share within each group, not at random, so
  # the panel spans the range of each class instead of sampling its middle.
  pick_ranks <- function(df, k = 8) {
    df <- df[order(df$dom_share), ]
    df[unique(round(seq(1, nrow(df), length.out = min(k, nrow(df))))), ]
  }
  vsel <- bind_rows(
    pick_ranks(rev[rev$crosses, ])  %>% mutate(class = "crosses"),
    pick_ranks(rev[!rev$crosses, ]) %>% mutate(class = "simple"))
  # The polygon the test reads is CLOSED, so the panel must draw the closing segment too: a
  # geom_path over the 100 samples alone leaves the loop open on screen and hides exactly the
  # crossings that involve the segment from the last sample back to the first.
  vtr <- readr::read_csv(ctf, show_col_types = FALSE) %>%
    mutate(stepIndex = strideIndex + stepInStride - 1) %>%
    inner_join(vsel %>% select(boutID, stepIndex, class, n_crossings, dom_share),
               by = c("boutID","stepIndex")) %>%
    group_by(boutID, stepIndex) %>%
    mutate(vfa_c = vfa_n - mean(vfa_n)) %>%
    arrange(pct, .by_group = TRUE) %>%
    slice(c(seq_len(dplyr::n()), 1L)) %>%
    mutate(panel = sprintf("%s %d, net %.0f%%", class, dplyr::first(n_crossings),
                           dplyr::first(dom_share))) %>%
    ungroup()
  pv <- ggplot(vtr, aes(vfa_c, vvert_n)) +
    geom_path(linewidth = 0.4, colour = "#3D4A5C") +
    geom_point(data = ~ dplyr::slice_head(dplyr::group_by(.x, panel), n = 1),
               size = 1.1, colour = "#C0398B") +
    facet_wrap(~ panel, ncol = 4) +
    coord_equal() + labs(x = LAB_FA_VEL, y = LAB_VERT_VEL) + theme_daley(base_size = 9) +
    theme(strip.text = element_text(size = 7))
  # One shared scale, because coord_equal() and free scales cannot both hold, and the shape of
  # the curve is what this panel is for. Canvas sized to the content so coord_equal() does not
  # pad the surplus width with white instead of enlarging the panels.
  ggsave("output_internal/hodograph_crossover_check.png", pv, width = 9, height = 8, dpi = 200)
  cat("  crossing-test validation panel: output_internal/hodograph_crossover_check.png\n")

  # per-gait rotation sense
  tab <- sa %>% group_by(gaitObjective) %>%
    summarise(n = n(), pct_CCW = round(100*mean(area > 0)), pct_CW = round(100*mean(area < 0)),
              median_area = round(median(area), 4), .groups = "drop")
  readr::write_csv(tab, "output/hodograph_rotation_by_gait.csv")
  cat("13_fig_rotation_sense: rotation sense by gait:\n"); print(as.data.frame(tab))

  # P(running-type = CW) vs dimensionless speed, with a bird random effect. Binomial GAM,
  # P(running-type) ~ s(u) + s(bird, bs="re"), with the crossover read as where the fitted
  # population probability crosses 0.5. The penalty is estimated from the data, so edf_u
  # reports how much curvature the speed term carries.
  have_mgcv <- requireNamespace("mgcv", quietly = TRUE)
  ucross <- NA_real_; gu <- NULL; note <- character(0)
  if (have_mgcv) {
    if (!exists("gam_reversal")) source("R/lib_gam_crossing.R")
    fitR <- gam_reversal(sa, sa$runType)      # shared model; see R/lib_gam_crossing.R
    ucross <- fitR$u_cross
    gu <- data.frame(meanSpeed_n = fitR$grid$u, p = fitR$grid$p_bouncing)
    note <- c("Fig 3 hodograph rotation-direction model.",
      "Signed area of the CoM velocity loop per step; sign = rotation sense (walking-type CCW > 0, running-type CW < 0).",
      "Binomial GAM P(running-type) ~ s(dimensionless speed u) + s(individual, random effect),",
      "over the individuals; fitted by the shared helper in R/lib_gam_crossing.R, which is also",
      "used for the energy-phase criterion in 12_hodograph_validation.R so the two crossover",
      "speeds are comparable and only one model definition exists.",
      sprintf("  population crossover speed (P=0.5) at u = %.2f (Froude %.2f); %d crossing(s) found.",
              ucross, ucross^2, fitR$n_crossings),
      sprintf(paste("  speed smooth effective degrees of freedom = %.2f of a maximum %g.",
                    "The smoothing penalty is estimated by REML, so the term shrinks towards a",
                    "straight line (edf 1) where the data support one; the fitted value is how",
                    "much curvature the reversal carries."),
              fitR$edf_u, fitR$edf_u_max),
      sprintf("By gait, %% walking-type (CCW): walk %d, grounded run %d, aerial run %d.",
              tab$pct_CCW[tab$gaitObjective=="walk"],
              tab$pct_CCW[tab$gaitObjective=="groundedRun"], tab$pct_CCW[tab$gaitObjective=="aerialRun"]),
      # Gait is classified from the aerial phase and the energy-phase criterion (02_clean.R),
      # neither of which uses the velocity loop, so these fractions are an observed
      # correspondence between two independent measures rather than a property of the labels.
      "Gait is classified from the aerial phase and KE-PE congruity, so rotation sense is measured independently of the labels and these fractions are an observed correspondence.")
    writeLines(note, "output_internal/hodograph_rotation_model_note.txt")
    readr::write_csv(data.frame(u = round(ucross, 4), Froude = round(ucross^2, 4),
                                n_crossings = fitR$n_crossings),
                     "output/rotation_reversal.csv")
    cat(sprintf("  P(running-type)=0.5 at dimensionless speed u=%.2f (Froude %.2f) [binomial GAM, %d crossing(s)].\n",
                ucross, ucross^2, fitR$n_crossings))
  }

  xlab_u <- LAB_U
  # Panel A: signed area vs speed, coloured by gait, zero line = walk/run boundary
  # The y title is set over two lines: at single-column width one line is longer than the
  # panel is tall and the axis label is clipped.
  pA <- ggplot(sa, aes(meanSpeed_n, area, colour = gaitObjective, shape = gaitObjective)) +
    geom_hline(yintercept = 0, linetype = 3, colour = "grey70") +
    geom_point(alpha = 0.6, size = 0.9) +
    scale_color_gait() + scale_shape_gait() +
    labs(x = xlab_u, y = expression(atop("Hodograph signed area",
                                         italic(A)/(italic(g) * italic(L)[0]) * "  (CCW +, CW -)"))) +
    theme_bio() + bio_drop_x() + bio_inset_legend(0.99, 0.02)
  # Panel B: P(running-type) vs speed with the binomial-GAM fit + empirical bins. With the
  # expanded sample the empirical proportions are sampled more finely (~100 steps
  # per bin, up to 25 bins) so the fit can be judged against a denser summary; each
  # point also carries a binomial 95% CI whisker.
  nb <- max(8, min(25, floor(nrow(sa)/100)))
  emp <- sa %>% mutate(ub = dplyr::ntile(meanSpeed_n, nb)) %>% group_by(ub) %>%
    summarise(u = mean(meanSpeed_n), p = mean(runType), k = sum(runType), n = dplyr::n(),
              lo = pmax(0, p - 1.96*sqrt(p*(1-p)/n)), hi = pmin(1, p + 1.96*sqrt(p*(1-p)/n)),
              .groups = "drop")
  pB <- ggplot(sa, aes(meanSpeed_n, runType))
  if (!is.null(gu)) {
    pB <- pB + geom_line(data = gu, aes(meanSpeed_n, p), colour = "grey20", linewidth = 1)
    if (is.finite(ucross)) pB <- pB + geom_vline(xintercept = ucross, linetype = 3, colour = "grey40")
  }
  pB <- pB +
    geom_errorbar(data = emp, aes(u, ymin = lo, ymax = hi), width = 0, colour = "grey55", inherit.aes = FALSE) +
    geom_point(data = emp, aes(u, p), size = 1.2, alpha = 0.85, inherit.aes = FALSE) +
    labs(x = xlab_u, y = "P(running-type rotation)") +
    theme_bio()

  # Both panels are drawn on the same speed axis, which panel B carries for the pair. The gait
  # key sits inside panel A, so the pair keeps the full single-column width for the data.
  fig <- (pA / pB) + plot_annotation(tag_levels = "A")
  bio_save("Fig3_RotationSense", fig, width = BIO_W1, height = 4.7)
  cat(sprintf("13_fig_rotation_sense: Fig3_RotationSense (%d steps).\n", nrow(sa)))
}
