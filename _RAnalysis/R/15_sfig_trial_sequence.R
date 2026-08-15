# 15_sfig_trial_sequence.R — Fig S4: every stride of one trial drawn as its own raw velocity
# loop, in travel order, so the rotation sense can be read stride by stride as the bird speeds up
# and slows down along the runway.
#
# Where Fig 5 averages many strides per gait after Hilbert registration, this draws the raw sampled
# trajectory of each individual stride with no registration and no averaging: fore-aft velocity
# fluctuation about that stride's mean on x, vertical velocity on y, both /sqrt(g L0), coloured by
# percent of the stride, with the stride touchdown a filled dot and the contralateral touchdown an
# open dot. Loops are laid out left to right in travel order, which is also speed order along each
# limb of the arc, so the sequence reads as the trial itself unfolds.
#
# Trial selection is by rule, not by hand, and the full ranking goes to
# output/trial_sequence_candidates.csv so the pick is auditable. A bout qualifies when it
#   (1) loses nothing anywhere in the pipeline at either grain (every detected step passes the QC
#       gate and enters the analysis sample, and every full step pair passes the stride QC gate),
#       and retains at least five full strides;
#   (2) uses at least two gaits;
#   (3) has a rise-then-fall speed arc: the peak speed falls strictly inside the trial, with at
#       least ARC_MIN_U of dimensionless speed gained before it and lost after it;
#   (4) changes gait on BOTH limbs of the arc: the gait at the first step and at the last step each
#       differ from the gait at the peak.
# Among qualifiers the pick is the most balanced arc (largest of min(rise, fall)), ties broken by
# mean vertical reconstruction drift.
#
# There is no steadiness gate: a trial with a speed arc is unsteady on its limbs by construction,
# so steadiness is annotated per stride instead. Reconstruction drift is judged against gait norms
# rather than across gaits, because drift scales strongly with gait
# (output/reconstruction_drift_by_gait.csv), so a walk-anchored trial cannot be compared with
# aerial-run trials on raw drift. A trailing odd step with no partner is drawn as a half loop and
# labelled from the step table.
#
# Outputs
#   output/SFig_TrialSequence.{pdf,png}
#   output/trial_sequence_candidates.csv   the ranked qualifying trials
#   output/trial_sequence_strides.csv      the drawn strides, with gait, sense, speed,
#                                          steadiness and the CW/CCW split of each loop's area
suppressPackageStartupMessages({library(dplyr); library(ggplot2)})
if (!exists("theme_daley")) source("R/lib_theme.R")
ctf <- "data/cycleTracesStep.csv"
if (!file.exists(ctf)) {
  cat("15_sfig_trial_sequence: cycleTracesStep.csv not found, skipping.\n")
} else {

ARC_MIN_U <- 0.3   # dimensionless speed that must be gained before the peak and lost after it

# The CW / CCW split of each loop's area comes from lib_hodograph.R, the one file that reads the
# closed velocity polygon, so the annotation on this figure and the crossover diagnostics in
# 13_fig_rotation_sense.R decompose the same quantity the classifier reads.
if (!exists("area_split")) source("R/lib_hodograph.R")

# `step` and `stride` are already the analysis sample (04_analysis_sample.R); step_raw and
# stride_raw give the pre-gate counts that the completeness rule compares against.
per_bout <- step %>%
  arrange(boutID, stepIndex) %>%
  group_by(boutID) %>%
  summarise(n_steps   = dplyr::n(),
            n_strides = { s0 <- 2 * floor((stepIndex - 1) / 2) + 1; sum(table(s0) == 2) },
            n_gaits  = n_distinct(as.character(gaitObjective)),
            gaits    = paste(unique(as.character(gaitObjective)), collapse = "+"),
            u_lo     = min(meanSpeed_n), u_hi = max(meanSpeed_n),
            peak_step  = which.max(meanSpeed_n),
            rise       = max(meanSpeed_n) - dplyr::first(meanSpeed_n),
            fall       = max(meanSpeed_n) - dplyr::last(meanSpeed_n),
            gait_start = dplyr::first(as.character(gaitObjective)),
            gait_peak  = as.character(gaitObjective)[which.max(meanSpeed_n)],
            gait_end   = dplyr::last(as.character(gaitObjective)),
            drift_mm = mean(driftRMS_vert_mm),
            we_resid = mean(abs(WE_relresid)),
            bird     = dplyr::first(as.character(bird)), .groups = "drop") %>%
  left_join(step_raw %>% count(boutID, name = "n_detected"), by = "boutID") %>%
  left_join(stride %>% count(boutID, name = "n_strides_qc"), by = "boutID") %>%
  mutate(n_strides_qc = coalesce(n_strides_qc, 0L),
         complete = n_steps == n_detected & n_strides_qc >= n_strides,
         arc_min  = pmin(rise, fall))

cand <- per_bout %>%
  filter(complete, n_strides >= 5, n_gaits >= 2,
         peak_step > 1, peak_step < n_steps,
         rise >= ARC_MIN_U, fall >= ARC_MIN_U,
         gait_start != gait_peak, gait_end != gait_peak) %>%
  arrange(desc(arc_min), drift_mm) %>%
  mutate(rank = dplyr::row_number(), .before = 1)
readr::write_csv(cand, "output/trial_sequence_candidates.csv")
stopifnot(nrow(cand) >= 1)
BOUT <- cand$boutID[1]

ct <- readr::read_csv(ctf, show_col_types = FALSE) %>% filter(boutID == BOUT)
stopifnot(nrow(ct) > 0)

loops <- ct %>%
  mutate(u = (stepInStride - 1) + pct / 100,          # 0..2 across the stride
         pctStride = 50 * u) %>%
  arrange(strideIndex, u) %>%
  group_by(strideIndex) %>%
  mutate(vfa_c = vfa_n - mean(vfa_n)) %>%
  ungroup() %>%
  mutate(seq = as.integer(factor(strideIndex)))       # 1..K in travel order

K  <- max(loops$seq)
dx <- 1.25 * max(tapply(loops$vfa_c, loops$seq, function(v) diff(range(v))))
loops <- loops %>% mutate(x = vfa_c + (seq - 1) * dx)
centers <- (seq_len(K) - 1) * dx

# Per-stride annotation. Full strides are labelled from the stride table; a trailing half stride
# has no stride row and is labelled from its own step's row instead.
half <- loops %>% group_by(strideIndex) %>%
  summarise(is_half = n_distinct(stepInStride) < 2, .groups = "drop")
step_fallback <- step %>% filter(boutID == BOUT) %>%
  transmute(strideIndex = stepIndex, gait_s = as.character(gaitObjective),
            hodo_s = hodoArea_n, u_s = meanSpeed_n, v_s = meanSpeed,
            acc_s = as.character(accClass))
asplit <- loops %>% arrange(strideIndex, u) %>% group_by(strideIndex) %>%
  summarise(v = list(area_split(vfa_c, vvert_n)), .groups = "drop") %>%
  mutate(pct_area_cw  = vapply(v, `[[`, numeric(1), "pct_area_cw"),
         pct_area_ccw = vapply(v, `[[`, numeric(1), "pct_area_ccw")) %>%
  select(-v)
# the two shares partition the loop's area
stopifnot(all(abs(asplit$pct_area_cw + asplit$pct_area_ccw - 100) < 1e-8))

info <- loops %>% distinct(strideIndex, seq) %>%
  left_join(asplit, by = "strideIndex") %>%
  left_join(stride %>% filter(boutID == BOUT) %>%
              select(strideIndex, gaitObjective, hodoArea_n, meanSpeed_n, meanSpeed, accClass),
            by = "strideIndex") %>%
  left_join(half, by = "strideIndex") %>%
  left_join(step_fallback, by = "strideIndex") %>%
  mutate(gait  = coalesce(as.character(gaitObjective), gait_s),
         hodo  = coalesce(hodoArea_n, hodo_s),
         u     = coalesce(meanSpeed_n, u_s),
         v_ms  = coalesce(meanSpeed, v_s),
         acc   = coalesce(as.character(accClass), acc_s),
         sense = ifelse(hodo > 0, "CCW", "CW"),
         lab1  = dplyr::recode(gait, walk = "walk", groundedRun = "gr. run",
                               aerialRun = "aer. run", .missing = "?"),
         extra = {e <- ifelse(acc == "steady", "",
                              dplyr::recode(acc, accelerating = "accel.", decelerating = "decel."))
                  e <- ifelse(is_half, trimws(paste(e, "half stride", sep = ", ")), e)
                  sub("^, ", "", e)},
         lab2  = sprintf("u = %.2f (%.2f m/s), %s\narea %.1f%% CW, %.1f%% CCW%s",
                         u, v_ms, sense, pct_area_cw, pct_area_ccw,
                         ifelse(extra == "", "", paste0("\n", extra))),
         col   = GAIT_COLORS[gait])
# The dominant share must reproduce the production sense, since both are the sign of the same sum.
stopifnot(all(ifelse(info$pct_area_ccw > 50, "CCW", "CW") == info$sense))
readr::write_csv(info %>% select(seq, strideIndex, gait, sense, u, v_ms, acc, is_half,
                                 pct_area_cw, pct_area_ccw),
                 "output/trial_sequence_strides.csv")

y_lo   <- min(loops$vvert_n); y_hi <- max(loops$vvert_n)
y_lab1 <- y_lo - 0.10; y_lab2 <- y_lo - 0.19
td  <- loops %>% filter(u == min(u))                          # stride touchdown
mid <- loops %>% filter(stepInStride == 2) %>% group_by(seq) %>%
  slice_min(u, n = 1, with_ties = FALSE) %>% ungroup()        # contralateral touchdown

p <- ggplot(loops, aes(x, vvert_n, group = seq)) +
  geom_hline(yintercept = 0, linetype = 3, colour = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = centers, linetype = 3, colour = "grey85", linewidth = 0.3) +
  geom_path(aes(colour = pctStride), linewidth = 0.9) +
  geom_point(data = td,  colour = "black", size = 1.8) +
  geom_point(data = mid, shape = 21, colour = "black", fill = "white", size = 1.8) +
  geom_text(data = info, aes(x = (seq - 1) * dx, y = y_lab1, label = lab1),
            colour = info$col, size = 2.6, fontface = "bold", inherit.aes = FALSE) +
  geom_text(data = info, aes(x = (seq - 1) * dx, y = y_lab2, label = lab2),
            colour = "grey40", size = 2.3, lineheight = 0.9, inherit.aes = FALSE) +
  scale_colour_gradientn(colours = c("#3D4A5C", "#1C9DA8", "#C0398B", "#D9A23B"),
                         name = "stride %", limits = c(0, 100)) +
  scale_x_continuous(breaks = centers, labels = seq_len(K)) +
  coord_equal() +
  labs(x = "stride within the trial (travel direction left to right)", y = LAB_VERT_VEL) +
  theme_daley()

# coord_equal(): size the canvas to the content so the panel is not padded with blank space
xr <- (K - 1) * dx + 1.6 * dx
yr <- (y_hi - y_lo) + 0.44                      # data height plus the three-line label band
S  <- min(4.5, 13 / xr)                         # inches per velocity unit, width-capped
W  <- xr * S + 1.7; H <- yr * S + 1.5
ggsave("output/SFig_TrialSequence.pdf", p, width = W, height = H, device = cairo_pdf)
ggsave("output/SFig_TrialSequence.png", p, width = W, height = H, dpi = 300, bg = "white")
cat(sprintf("15_sfig_trial_sequence: %d bouts pass the arc rule; chose %s (%s), %d strides.\n",
            nrow(cand), BOUT, cand$bird[1], K))
print(as.data.frame(info %>% select(seq, gait, sense, u, acc, is_half)), digits = 3)
}
