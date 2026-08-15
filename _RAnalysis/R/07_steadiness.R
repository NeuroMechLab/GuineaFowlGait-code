# 07_steadiness.R (supplement) — steadiness classification, its speed coverage, and the evidence
# that the criterion is speed-neutral. Criterion (accClass): net CoM mechanical energy change per
# unit distance as a fraction of body weight, |fracEG| = |dE_CoM/(m g L)| <= STEADY_GRADE, an
# effective fore-aft grade. Axis: dimensionless speed u = v/sqrt(gL0).
#
# The manuscript reports only the criterion it uses. The comparison against the two alternatives
# defined in 02_clean.R (accClass_en, |dE_CoM/(m v^2)| <= STEADY_FRAC_E; accClass_fa,
# |a_fa T/v| <= STEADY_FRAC) is still computed below and written to
# output_internal/steadiness_criteria_comparison.csv, because it is the internal record behind the choice
# of criterion, but it is NOT plotted and no manuscript
# number is drawn from it. Speed-neutrality is shown positively instead, as the exclusion rate of
# the criterion actually used, flat across the speed range (panel C).
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2); library(patchwork)})
if (!dir.exists("output")) dir.create("output")
xlab_u <- LAB_U

# ---- counts + criterion comparison ---------------------------------------
stride_counts <- stride %>% count(accClass, name = "strides") %>%
  mutate(pct = round(100*strides/sum(strides), 1))
step_counts   <- step %>% count(accClass, name = "steps") %>%
  mutate(pct = round(100*steps/sum(steps), 1))
readr::write_csv(stride_counts, "output/steadiness_stride_counts.csv")
readr::write_csv(step_counts,   "output/steadiness_step_counts.csv")
cmp <- stride %>% summarise(
  grade_steady_pct      = round(100*mean(accClass_eg == "steady")),
  energyfrac_steady_pct = round(100*mean(accClass_en == "steady")),
  foreaft_steady_pct    = round(100*mean(accClass_fa == "steady")),
  grade_vs_energyfrac_agreement_pct = round(100*mean(accClass_eg == accClass_en)))
readr::write_csv(cmp, "output_internal/steadiness_criteria_comparison.csv")
cat(sprintf("Steadiness: PRIMARY = energy grade |dE_CoM/(m g L)|<=%.2f; previous energy-fraction |fracDE|<=%.2f; fore-aft |fracDV|<=%.2f.\n",
            STEADY_GRADE, STEADY_FRAC_E, STEADY_FRAC))
cat("Per-stride (primary/energy grade):\n"); print(stride_counts)
cat(sprintf("Criteria comparison (strides): grade %d%%, previous energy-fraction %d%%, fore-aft %d%% steady; grade-vs-previous agreement %d%%.\n",
            cmp$grade_steady_pct, cmp$energyfrac_steady_pct, cmp$foreaft_steady_pct,
            cmp$grade_vs_energyfrac_agreement_pct))

# per-gait steady/accel/decel strides, so the steady sample behind Fig 5 is explicit per gait
if ("gaitObjective" %in% names(stride)) {
  by_gait <- stride %>% filter(!is.na(gaitObjective)) %>%
    count(gaitObjective, accClass) %>%
    tidyr::pivot_wider(names_from = accClass, values_from = n, values_fill = 0)
  readr::write_csv(by_gait, "output/steadiness_by_gait.csv")
  cat("Steady/accel/decel strides by objective gait (primary criterion):\n")
  print(as.data.frame(by_gait))
}

# ---- speed coverage: steady subset vs full ------------------------------
# Reported at BOTH grains. The criterion is evaluated per stride, so the figure is drawn at
# stride grain (see the note above panel B); the step-grain row is kept because the steady
# subset used by Table 1 and Figs 4, 6 and 7 is a set of steps.
cov_of <- function(d, grain) {
  sd_ <- d %>% filter(accClass == "steady")
  tibble::tibble(grain = grain, subset = c("full","steady"), n = c(nrow(d), nrow(sd_)),
    u_min = c(min(d$meanSpeed_n), min(sd_$meanSpeed_n)),
    u_p05 = c(quantile(d$meanSpeed_n, .05), quantile(sd_$meanSpeed_n, .05)),
    u_med = c(median(d$meanSpeed_n), median(sd_$meanSpeed_n)),
    u_max = c(max(d$meanSpeed_n), max(sd_$meanSpeed_n)))
}
cov <- bind_rows(cov_of(step, "step"), cov_of(stride, "stride"))
readr::write_csv(cov, "output/steadiness_speed_coverage.csv")
steady <- stride %>% filter(accClass == "steady")

# ---- figure --------------------------------------------------------------
# steadiness palette, deliberately distinct from the gait palette (navy/teal/magenta)
pal <- c(decelerating = "#8E63A6", steady = "#555555", accelerating = "#E69F00")
# Every panel is at STRIDE grain, because the criterion is evaluated per stride. Plotting strides
# shows the criterion the paper actually uses.
pA <- ggplot() +
  geom_density(data = stride, aes(meanSpeed_n, after_stat(count), fill = "full"),   alpha = 0.35, colour = NA) +
  geom_density(data = steady, aes(meanSpeed_n, after_stat(count), fill = "steady"), alpha = 0.55, colour = NA) +
  scale_fill_manual(values = c(full = "#B0B7C0", steady = "#3D4A5C"),
                    labels = c("all strides","steady only"), name = NULL) +
  labs(x = xlab_u, y = "stride count (density-scaled)") + theme_daley()

# B. the primary criterion itself: energy grade vs speed with the +/- band
pB <- ggplot(stride, aes(meanSpeed_n, 100*fracEG, colour = accClass)) +
  geom_hline(yintercept = c(-100*STEADY_GRADE, 100*STEADY_GRADE), linetype = 3, colour = "grey50") +
  geom_point(alpha = 0.5, size = 1.3) +
  scale_colour_manual(values = pal, name = NULL) +
  labs(x = xlab_u, y = "net CoM energy change per distance / body weight (%)") +
  theme_daley()

# C. speed-neutrality of the criterion, stated positively: the fraction of steps excluded as
# unsteady, in equal-count speed bins. A flat line is the property the Methods claim, that the
# criterion does not preferentially exclude slow or fast steps.
ubrk <- quantile(stride$meanSpeed_n, seq(0, 1, length.out = 7), na.rm = TRUE)
excl <- stride %>%
  mutate(ubin = cut(meanSpeed_n, ubrk, include.lowest = TRUE), ex = accClass != "steady") %>%
  group_by(ubin) %>%
  summarise(u = median(meanSpeed_n), rate = 100*mean(ex), n = dplyr::n(), .groups = "drop")
readr::write_csv(excl %>% select(u, exclusion_pct = rate, n_strides = n),
                 "output/steadiness_exclusion_by_speed.csv")
pC <- ggplot(excl, aes(u, rate)) +
  geom_line(linewidth = 1, colour = "#1C9DA8") + geom_point(size = 2, colour = "#1C9DA8") +
  ylim(0, 100) +
  labs(x = xlab_u, y = "strides excluded as unsteady (%)") +
  theme_daley()

fig <- (pA / pB / pC) + plot_layout(guides = "collect") +
  patchwork::plot_annotation(tag_levels = "A")
ggsave("output/SFig_Steadiness.pdf", fig, width = 7, height = 11, device = cairo_pdf)
ggsave("output/SFig_Steadiness.png", fig, width = 7, height = 11, dpi = 300)
cat("07_steadiness: wrote counts, criterion comparison, coverage and SFig_Steadiness.\n")
