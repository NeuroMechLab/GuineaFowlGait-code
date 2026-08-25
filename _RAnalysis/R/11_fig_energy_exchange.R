# 11_fig_energy_exchange.R — Fig 4: how the traditional and proposed CoM-dynamics
# metrics map onto the classified gaits, vs dimensionless speed (u = v/sqrt(gL0)),
# gait-coded, steady STEPS only. Three descriptors, none used as a classifier:
# pendular energy recovery (Cavagna et al. 1977), KE-PE congruity (Ahn, Furrow &
# Biewener 2004), and the collision angle Phi (Lee et al. 2011, 2013; Phi is
# equivalent to the dimensionless mechanical cost of transport). The point is the
# continuous variation within as well as between gaits.
suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(patchwork)})
if (!exists("theme_daley")) source("R/lib_theme.R")
gcol <- "gaitObjective"
xlab <- LAB_U
d <- step %>% filter(accClass == "steady", is.finite(recovery), is.finite(congruity))

# Collision angle vs mechanical cost of transport correspondence (all classified steps).
# Reported in the paper as evidence that the force-and-velocity-weighted collision angle
# is equivalent to the dimensionless mechanical cost of transport (Lee et al. 2013).
cc <- step %>% filter(gaitObjective %in% GAIT_LEVELS, is.finite(collisionAngle), is.finite(CoTmech))
r_ccot <- cor(cc$collisionAngle, cc$CoTmech)
# Cost of transport tracks the collision angle and not pendular recovery, so both correlations are
# emitted. The recovery pair is the same quantity Table S2 prints from descriptor_correlations.csv;
# it is computed here for the within-bird companion below.
cr <- step %>% filter(gaitObjective %in% GAIT_LEVELS, is.finite(CoTmech), is.finite(recovery))
r_crec <- cor(cr$CoTmech, cr$recovery)
# Pooled correlations treat every step as independent, which they are not, since each bird
# contributes many. For a
# correlation carried by differences BETWEEN birds rather than by any relationship within one, the
# pooled value can take the opposite sign to the typical within-bird value, as it does for cost of
# transport against recovery. The within-bird distribution is therefore emitted for every pooled
# value, to output_internal/ since the paper reports the pooled values only.
within_bird <- function(dat, x, y, lab, min_n = 30) {
  w <- dat %>% group_by(subjectID) %>% filter(dplyr::n() >= min_n) %>%
    summarise(r = cor(.data[[x]], .data[[y]]), n = dplyr::n(), .groups = "drop")
  data.frame(metric = lab, n_birds = nrow(w),
             r_within_median = round(median(w$r), 3),
             r_within_min = round(min(w$r), 3), r_within_max = round(max(w$r), 3))
}
readr::write_csv(data.frame(
    metric = c("collision_cot_pearson_r", "cot_recovery_pearson_r"),
    value  = c(round(r_ccot, 3), round(r_crec, 3)),
    n      = c(nrow(cc), nrow(cr))),
  "output/energy_exchange_stats.csv")
readr::write_csv(bind_rows(
    within_bird(cc, "collisionAngle", "CoTmech", "collision_cot_pearson_r"),
    within_bird(cr, "CoTmech", "recovery", "cot_recovery_pearson_r")),
  "output_internal/energy_exchange_within_bird.csv")
cat("11_fig_energy_exchange: within-bird correlation spread:\n")
print(as.data.frame(bind_rows(
  within_bird(cc, "collisionAngle", "CoTmech", "collision_cot_pearson_r"),
  within_bird(cr, "CoTmech", "recovery", "cot_recovery_pearson_r"))))
cat(sprintf("11_fig_energy_exchange: collision-angle ~ CoT Pearson r = %.3f (n = %d steps); CoT ~ recovery r = %.3f (n = %d)\n",
            r_ccot, nrow(cc), r_crec, nrow(cr)))

# Trend fit: a penalized cubic-regression spline (mgcv GAM, y ~ s(u)) fit on
# EQUAL-WIDTH speed-bin means rather than a loess through the raw points. Loess (and
# a GAM on the raw points) is dominated by the mid-speed range where sample density
# is highest, so it fits the sparse slow and fast edges poorly. Binning to equal
# speed intervals first gives every speed region equal weight, so the trend is not
# pulled by the dense middle; bin means (with 95% CI whiskers) are drawn so the
# reader sees the unbiased summary the curve is fit to.
# The three panels share one speed axis, carried by the lowest of them, so all three are drawn
# on the same x range. ymax caps a panel's view: coord_cartesian sets the view, so a point above
# the cap is off the panel and still enters the bin means and the spline.
XLIM <- range(d$meanSpeed_n, na.rm = TRUE)
metric_panel <- function(y, ylab, ymax = NA) {
  dd <- d %>% filter(is.finite(.data[[y]]))
  rng <- range(dd$meanSpeed_n, na.rm = TRUE)
  brk <- seq(rng[1], rng[2], length.out = 19)                 # 18 equal-width speed bins
  bn <- dd %>% mutate(sb = cut(meanSpeed_n, brk, include.lowest = TRUE)) %>%
    group_by(sb) %>% summarise(u = mean(meanSpeed_n), m = mean(.data[[y]]),
      se = sd(.data[[y]]) / sqrt(dplyr::n()), n = dplyr::n(), .groups = "drop") %>%
    filter(n >= 3)
  ggplot(dd, aes(meanSpeed_n, .data[[y]])) +
    geom_point(aes(colour = .data[[gcol]], shape = .data[[gcol]]), alpha = 0.55, size = 0.9) +
    geom_errorbar(data = bn, aes(u, ymin = m - 1.96*se, ymax = m + 1.96*se),
                  width = 0, colour = "grey35", inherit.aes = FALSE, linewidth = 0.3) +
    geom_point(data = bn, aes(u, m), colour = "grey15", size = 0.9, inherit.aes = FALSE) +
    geom_smooth(data = bn, aes(u, m), method = "gam",
                formula = y ~ s(x, bs = "cs", k = 6), se = TRUE, colour = "grey15", linewidth = 0.5) +
    scale_color_gait() + scale_shape_gait() +
    coord_cartesian(xlim = XLIM, ylim = if (is.na(ymax)) NULL else c(NA, ymax)) +
    labs(x = xlab, y = ylab) + theme_bio()
}
p_rec <- metric_panel("recovery",  "Pendular energy recovery (%)")
p_con <- metric_panel("congruity", "KE-PE congruity (%)")
# Collision angle is plotted in DEGREES, the form reported in the tables and the text. The
# theoretical equivalence with the dimensionless mechanical cost of transport holds for the angle
# in radians, so the axis is not labelled with that equivalence; it is stated in the caption
# and quantified by the correlation above, which is computed on the radian form.
ycol  <- if ("collisionAngle_deg" %in% names(d)) "collisionAngle_deg" else "CoTmech"
# Capped so the panel resolves the band the steps occupy; the steps above the cap are counted
# into the artifact the legend quotes.
COLL_CAP <- 25
p_col <- metric_panel(ycol, "Collision angle (deg)", ymax = COLL_CAP)
readr::write_csv(
  data.frame(panel = "Fig 4C", measure = ycol, axis_cap_deg = COLL_CAP,
             steps_drawn = sum(is.finite(d[[ycol]])),
             steps_above_cap = sum(d[[ycol]] > COLL_CAP, na.rm = TRUE),
             max_deg = round(max(d[[ycol]], na.rm = TRUE), 1)),
  "output/fig4_axis_clipping.csv")

# The gait key sits inside panel B, in the corner the congruity cloud leaves empty.
NOKEY <- theme(legend.position = "none")
fig <- (p_rec + bio_drop_x() + NOKEY) /
       (p_con + bio_drop_x() + bio_inset_legend(0.99, 0.02)) /
       (p_col + NOKEY) +
  plot_annotation(tag_levels = "A")
bio_save("Fig4_MetricMapping_Steady", fig, width = BIO_W1, height = 6.4)
cat(sprintf("11_fig_energy_exchange: Fig4_MetricMapping_Steady (%d steady steps).\n", nrow(d)))

# ---- how closely does the weighted collision angle track the mechanical cost of transport? ----
# Lee et al. (2011, equations 6 to 12) derive the link as a TWO-STAGE APPROXIMATION, not an
# identity:
#   (1) the small-angle substitution sin(phi) ~ phi, which they state holds for phi below about
#       0.3 rad, taking the weighted collision angle Phi to the mechanical cost of motion; and
#   (2) mean|V| ~ mean forward velocity (their eq 9) and mean|F| ~ body weight (their eq 10),
#       which hold for small vertical and lateral oscillations and small fore-aft and lateral
#       forces, taking the cost of motion to the cost of transport.
# Algebraically, with sin(phi) = |F.V|/(|F||V|) by definition, Phi/CoT = (BW vbar)/mean(|F||V|)
# under (1) alone, so stage (2) is exactly the statement that this ratio is one. That is emitted
# here per gait, so the paper can say how well the approximation holds in these data rather than
# asserting an equivalence. The step-to-step spread of the ratio is what limits the correlation.
if (file.exists("data/cycleTracesStep.csv")) {
  ce <- readr::read_csv("data/cycleTracesStep.csv", show_col_types = FALSE) %>%
    mutate(stepIndex = strideIndex + stepInStride - 1) %>%
    inner_join(step %>% select(boutID, stepIndex, gaitObjective, collisionAngle, CoTmech,
                               meanSpeed_n),
               by = c("boutID", "stepIndex")) %>%
    group_by(boutID, stepIndex, gaitObjective, collisionAngle, CoTmech, meanSpeed_n) %>%
    summarise(meanFmag_BW = mean(sqrt(Fz_BW^2 + Ffa_BW^2)),
              meanVmag_n  = mean(sqrt(vfa_n^2 + vvert_n^2)),
              meanFV      = mean(sqrt(Fz_BW^2 + Ffa_BW^2) * sqrt(vfa_n^2 + vvert_n^2)),
              .groups = "drop") %>%
    mutate(vbar = abs(meanSpeed_n),
           eq10_meanF_over_BW   = meanFmag_BW,
           eq09_meanV_over_vbar = meanVmag_n / vbar,
           product_over_BWvbar  = meanFV / vbar,
           ratio_Phi_over_CoT   = collisionAngle / CoTmech)
  cc_tab <- ce %>% group_by(gait = gaitObjective) %>%
    summarise(n_steps = dplyr::n(),
              median_Phi_rad         = round(median(collisionAngle), 3),
              pct_Phi_above_0.3rad   = round(100 * mean(collisionAngle > 0.3), 1),
              median_eq10_meanF_BW   = round(median(eq10_meanF_over_BW), 3),
              median_eq09_meanV_vbar = round(median(eq09_meanV_over_vbar), 3),
              median_ratio_Phi_CoT   = round(median(ratio_Phi_over_CoT), 3),
              ratio_q25 = round(quantile(ratio_Phi_over_CoT, .25), 3),
              ratio_q75 = round(quantile(ratio_Phi_over_CoT, .75), 3),
              pearson_r_Phi_CoT      = round(cor(collisionAngle, CoTmech), 3),
              .groups = "drop")
  overall <- ce %>% summarise(gait = "all steps", n_steps = dplyr::n(),
              median_Phi_rad = round(median(collisionAngle), 3),
              pct_Phi_above_0.3rad = round(100 * mean(collisionAngle > 0.3), 1),
              median_eq10_meanF_BW = round(median(eq10_meanF_over_BW), 3),
              median_eq09_meanV_vbar = round(median(eq09_meanV_over_vbar), 3),
              median_ratio_Phi_CoT = round(median(ratio_Phi_over_CoT), 3),
              ratio_q25 = round(quantile(ratio_Phi_over_CoT, .25), 3),
              ratio_q75 = round(quantile(ratio_Phi_over_CoT, .75), 3),
              pearson_r_Phi_CoT = round(cor(collisionAngle, CoTmech), 3))
  readr::write_csv(bind_rows(overall, cc_tab), "output/collision_cot_relationship.csv")
  cat("11_fig_energy_exchange: collision angle vs cost of transport, Lee et al. approximation:\n")
  print(as.data.frame(bind_rows(overall, cc_tab)))
}

# ---- is the enclosed hodograph area a measure of the mechanical cost of transport? ----------
# The collision framing links the step-to-step velocity redirection to the work of transport, so
# the magnitude of the loop area is a candidate cost measure alongside the collision angle. Test
# it directly: |A|/(gL0) against the dimensionless mechanical cost of transport, over the
# analysis sample and within each gait, with the collision angle as the reference comparison.
# Spearman is emitted alongside Pearson because the area magnitude is right-skewed.
ha <- step %>% filter(gaitObjective %in% GAIT_LEVELS, is.finite(hodoArea_n), is.finite(CoTmech),
                      is.finite(collisionAngle))
area_cot <- function(dat, lab) data.frame(
  group = lab, n_steps = nrow(dat),
  r_absArea_CoT      = round(cor(abs(dat$hodoArea_n), dat$CoTmech), 3),
  rho_absArea_CoT    = round(cor(abs(dat$hodoArea_n), dat$CoTmech, method = "spearman"), 3),
  r_signedArea_CoT   = round(cor(dat$hodoArea_n, dat$CoTmech), 3),
  r_absArea_collAng  = round(cor(abs(dat$hodoArea_n), dat$collisionAngle), 3),
  r_collAng_CoT      = round(cor(dat$collisionAngle, dat$CoTmech), 3),
  median_absArea     = round(median(abs(dat$hodoArea_n)), 4),
  median_CoT         = round(median(dat$CoTmech), 3))
area_cot_tbl <- bind_rows(
  area_cot(ha, "all steps"),
  do.call(rbind, lapply(GAIT_LEVELS, function(g) area_cot(ha[ha$gaitObjective == g, ], g))))
readr::write_csv(area_cot_tbl, "output/hodoarea_cot_relationship.csv")
cat("11_fig_energy_exchange: hodograph |area| vs mechanical cost of transport:\n")
print(as.data.frame(area_cot_tbl))
