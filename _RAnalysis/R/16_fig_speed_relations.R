# 16_fig_speed_relations.R — Fig 6: four-row panel vs DIMENSIONLESS SPEED
# u = v/sqrt(gL0) (steady STEPS only): duty factor, normalized step length,
# normalized step frequency, and the observed distribution of steps over speed in each of the
# three gaits.
#
# Panel D shows the speed distributions directly, as histograms of the observed steps per gait.
# Gait is assigned from two qualitative features of the same steps, so the counts already carry
# the result, which is how far the gaits overlap in speed.
suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(patchwork)})
if (!exists("theme_daley")) source("R/lib_theme.R")
gcol <- "gaitObjective"
d <- step %>% filter(accClass == "steady")

xlab_u <- LAB_U
# Trend fit: a penalized cubic-regression spline (GAM) fit on EQUAL-WIDTH speed-bin means, not
# a loess through the raw points, so the sparse slow and fast edges are not overweighted by the
# dense mid-speed range (same rationale and method as Fig 3).
scat <- function(y, ylab, href = NA, xlab = NULL) {
  dd <- d %>% filter(is.finite(.data[[y]]))
  brk <- seq(min(dd$meanSpeed_n), max(dd$meanSpeed_n), length.out = 19)
  bn <- dd %>% mutate(sb = cut(meanSpeed_n, brk, include.lowest = TRUE)) %>%
    group_by(sb) %>% summarise(u = mean(meanSpeed_n), m = mean(.data[[y]]),
      n = dplyr::n(), .groups = "drop") %>% filter(n >= 3)
  p <- ggplot(dd, aes(meanSpeed_n, .data[[y]])) +
    geom_point(aes(colour = .data[[gcol]], shape = .data[[gcol]]), alpha = 0.55, size = 0.8) +
    geom_smooth(data = bn, aes(u, m), method = "gam", formula = y ~ s(x, bs = "cs", k = 6),
                se = TRUE, colour = "grey20", linewidth = 0.5) +
    scale_color_gait() + scale_shape_gait() + labs(x = xlab, y = ylab) +
    coord_cartesian(xlim = XLIM) + theme_bio()
  if (is.finite(href)) p <- p + geom_hline(yintercept = href, linetype = 3, colour = "grey55")
  p
}
# The four panels share one speed axis, carried by panel D, so all four are drawn on the same
# x range.
XLIM <- range(d$meanSpeed_n, na.rm = TRUE)
pA <- scat("dutyFactor",  "Duty factor", href = 0.5)
pB <- scat("stepLength_n", LAB_STEPLEN)
pC <- scat("stepFreq_n", LAB_STEPFREQ)

# Panel D: observed distribution of steady steps over dimensionless speed, one histogram per
# gait. Bin width is common to the three gaits so the counts are directly comparable, and each
# gait is drawn with its own outline THICKNESS as well as its colour, so the panel is readable
# without colour and without a dash pattern. Counts, not densities, because the relative size of the three samples is part
# of what the panel shows.
ds <- d %>% filter(as.character(gaitObjective) %in% GAIT_LEVELS) %>%
  mutate(gait = factor(as.character(gaitObjective), levels = GAIT_LEVELS))
BINW <- 0.1
hist_tbl <- ds %>%
  mutate(bin = BINW * (floor(meanSpeed_n / BINW) + 0.5)) %>%
  count(gait, u_bin_center = bin, name = "steps")
readr::write_csv(hist_tbl, "output/gait_speed_histogram_steady.csv")
ovl <- ds %>% group_by(gait) %>%
  summarise(n = dplyr::n(), u_min = min(meanSpeed_n), u_max = max(meanSpeed_n),
            u_q25 = quantile(meanSpeed_n, 0.25), u_med = median(meanSpeed_n),
            u_q75 = quantile(meanSpeed_n, 0.75), .groups = "drop") %>%
  mutate(across(where(is.numeric) & !n, ~round(.x, 2)))
readr::write_csv(ovl, "output/gait_speed_overlap_steady.csv")

pD <- ggplot(ds, aes(meanSpeed_n)) +
  geom_histogram(aes(fill = gait, colour = gait, linewidth = gait),
                 binwidth = BINW, position = "identity", alpha = 0.35) +
  scale_fill_gait() + scale_color_gait() +
  scale_linewidth_manual(values = GAIT_LWD, labels = GAIT_LABELS,
                         breaks = GAIT_LEVELS, name = "Gait") +
  labs(x = xlab_u, y = "Steps per speed bin") +
  coord_cartesian(xlim = XLIM) +                       # share the speed axis with A-C
  theme_bio() +
  guides(colour = "none", linewidth = "none", fill = "none")  # one Gait legend, from panels A-C
cat("Speed overlap between gaits, steady steps:\n"); print(as.data.frame(ovl))

# The gait key sits inside panel B, in the corner the step-length trend leaves empty.
NOKEY <- theme(legend.position = "none")
fig <- (pA + bio_drop_x() + NOKEY) /
       (pB + bio_drop_x() + bio_inset_legend(0.99, 0.02)) /
       (pC + bio_drop_x() + NOKEY) /
       (pD + NOKEY) +
  plot_annotation(tag_levels = "A")
bio_save("Fig6_SpeedRelations_Steady", fig, width = BIO_W1, height = 7.6)
cat(sprintf("16_fig_speed_relations: Fig6_SpeedRelations_Steady (%d steady steps).\n", nrow(d)))
