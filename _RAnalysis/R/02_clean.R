# 02_clean.R — apply the QC gate, factors, gait ordering, and steadiness labels.
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(readr)})

# QC-retention gate. The gate (per-observation work-energy residual and duration-normalized
# vertical CoM drift) and the per-step energy steadiness are computed ONCE in MATLAB
# (GaitSelMulti_QCGate) and written to data/step_qcpass.csv, data/stride_qcpass.csv and
# data/step_steadiness.csv. That is the single source of truth, so the Dryad package, the mean
# traces, and these figures share one definition. Here we read the flags and keep qcPass == 1.
step_qc   <- readr::read_csv("data/step_qcpass.csv",   show_col_types = FALSE)
stride_qc <- readr::read_csv("data/stride_qcpass.csv", show_col_types = FALSE)

clean_common <- function(df, qc, keycol) {
  df %>%
    left_join(qc, by = c("boutID", keycol)) %>%
    filter(qcPass == 1) %>%
    select(-qcPass) %>%
    mutate(bird = factor(bird),
           gaitHand = factor(gaitLabelHand, levels = c("walk","grounded","aerial")),
           accelDir = factor(ifelse(accelForeAft >= 0, "accel", "decel"),
                             levels = c("decel","accel")))
}

stride <- clean_common(stride_raw, stride_qc, "strideIndex")
step   <- clean_common(step_raw, step_qc, "stepIndex")
cat(sprintf("After QC (MATLAB gate): %d strides, %d steps retained.\n", nrow(stride), nrow(step)))

# The MATLAB export's gait label is a function of the hodograph signed area (see the
# classification block below). It is renamed here to gaitHodo so that the dependency of any
# result on the hodograph stays legible in the code, and gaitObjective is rebuilt from
# features the hodograph does not enter. gaitHodo is used only by 12_hodograph_validation.R.
rename_hodo_label <- function(df) {
  if ("gaitObjective" %in% names(df))
    df <- df %>% rename(gaitHodo = gaitObjective) %>%
      mutate(gaitHodo = factor(gaitHodo, levels = GAIT_LEVELS))
  if ("mechClass" %in% names(df)) df <- df %>% rename(mechHodo = mechClass)
  df
}
stride <- rename_hodo_label(stride); step <- rename_hodo_label(step)

# Steadiness classification. The PRIMARY criterion is the net change in total CoM
# mechanical energy per unit distance, as a fraction of body weight,
#   fracEG = dE_CoM / (m g L),   L = strideLength (stride) or stepLength (step),
# classified steady when |fracEG| <= STEADY_GRADE. This is the mean net fore-aft force
# over the stride relative to body weight, equivalently the effective fore-aft grade the
# CoM energy climbs or descends (STEADY_GRADE = 0.05 is a 5% grade, ~2.9 deg). It is
# dimensionless, size-neutral and speed-neutral: both the numerator and the distance scale grow
# with speed. The exclusion rate it produces is nonetheless not flat: it is emitted to
# steadiness_exclusion_by_speed.csv and rises in the fastest bin, which 07_steadiness.R panel C
# shows and the Methods now states.
# The energy-fraction form |dE_CoM/(m v^2)| is kept for comparison as accClass_en.
# The fore-aft speed-change criterion (fracDV, Birn-Jeffery & Daley 2012; Birn-Jeffery et al.
# 2014) is also kept for comparison. fracEG is evaluated per stride (its touchdown-cut-phase
# component cancels over a stride) and propagated to the two steps of the stride; steps without
# a matched stride fall back to their own fracEG. Must match GaitSelMulti_QCGate.m.
STEADY_GRADE  <- 0.05   # primary: |dE_CoM/(m g L)| effective fore-aft grade
STEADY_FRAC   <- 0.10   # comparison: fore-aft |a_fa*T/v|; +/-1 s.d. of level trials (~10% speed)
STEADY_FRAC_E <- 0.10   # comparison: previous energy fraction |dE_CoM/(m v^2)|
G_ACC <- 9.81
classify_acc <- function(x, thr) factor(ifelse(x > thr, "accelerating",
                                        ifelse(x < -thr, "decelerating", "steady")),
                                        levels = c("decelerating","steady","accelerating"))
stride <- stride %>% mutate(
  fracEG      = dE_CoM/(mass_kg*G_ACC*strideLength),       # primary metric (signed)
  fracDV      = accelForeAft*stridePeriod/meanSpeed,
  accClass_eg = classify_acc(fracEG, STEADY_GRADE),
  accClass_fa = classify_acc(fracDV, STEADY_FRAC),
  accClass_en = classify_acc(fracDE, STEADY_FRAC_E),
  accClass    = accClass_eg)                               # primary = energy grade per stride
step <- step %>% mutate(
  fracEG        = dE_CoM/(mass_kg*G_ACC*stepLength),
  fracDV        = accelForeAft*stepPeriod/meanSpeed,
  accClass_eg   = classify_acc(fracEG, STEADY_GRADE),
  accClass_fa   = classify_acc(fracDV, STEADY_FRAC),
  accClass_en   = classify_acc(fracDE, STEADY_FRAC_E),
  accClass_step = accClass_eg,
  s0            = 2*floor((stepIndex-1)/2) + 1)
step <- step %>%
  left_join(stride %>% select(boutID, strideIndex, accClass_par = accClass_eg),
            by = c("boutID", "s0" = "strideIndex")) %>%
  mutate(accClass = dplyr::coalesce(accClass_par, accClass_step)) %>%   # inherit stride steadiness
  select(-accClass_par)
# ---- Gait classification from two qualitative features, neither of them the hodograph ----
# Feature 1, the aerial phase: the summed vertical force reaches zero (hasFlight). This is
# the clean form of the duty-factor walk-run criterion.
#
# Feature 2, pendular against bouncing energy exchange: the CLASSIC Cavagna criterion in its
# per-cycle form. KE and PE change in the same direction over a fraction of the cycle
# (congruity, Ahn et al. 2004); out of phase for most of the cycle is pendular and in phase
# for most of it is bouncing, so the boundary is 50%, the in-phase against out-of-phase
# divide, rather than a tuned level. A few steps sit exactly at 50% and are classified bouncing
# by the strict inequality; the count is emitted below so the tie rule is visible.
#
# Congruity is computed from the CoM kinetic and potential energies and never from the velocity
# loop, so the hodograph is independent of this classification and can be tested against it in
# 12_hodograph_validation.R. The hodograph-based label travels alongside as gaitHodo
# (GaitSel_PerStepStrideMeasures/classify), because the paper reports both.
CONGRUITY_SPLIT <- 50   # percent of the cycle with KE and PE changing in the same direction

gait4_of <- function(hasFlight, pendular) factor(dplyr::case_when(
  hasFlight == 0 &  pendular ~ "walkGrounded",
  hasFlight == 0 & !pendular ~ "groundedRun",
  hasFlight == 1 & !pendular ~ "aerialRun",
  hasFlight == 1 &  pendular ~ "pendularRun"), levels = GAIT4_LEVELS)

add_gait <- function(df) df %>% mutate(
  pendular  = congruity < CONGRUITY_SPLIT,
  mechClass = factor(ifelse(pendular, "pendular", "bouncing"),
                     levels = c("pendular", "bouncing")),
  # three-gait form used by every figure: every step with a flight phase is aerial running
  gaitObjective = factor(dplyr::case_when(
    hasFlight == 0 &  pendular ~ "walk",
    hasFlight == 0 & !pendular ~ "groundedRun",
    hasFlight == 1             ~ "aerialRun"), levels = GAIT_LEVELS),
  # four-cell form used by the summary tables, which splits the aerial-pendular cell out
  gait4     = gait4_of(hasFlight, pendular),
  # the hodograph-based four-cell form, for the validation script only
  gait4Hodo = gait4_of(hasFlight, hodoArea >= 0))
step <- add_gait(step); stride <- add_gait(stride)

cat("Four-cell gait counts (steps), classic energy-phase criterion:\n")
print(table(step$gait4, useNA = "ifany"))
cat(sprintf("Steps with congruity exactly %g%% (tie, classified bouncing): %d\n",
            CONGRUITY_SPLIT, sum(step$congruity == CONGRUITY_SPLIT, na.rm = TRUE)))
if ("gaitHodo" %in% names(step)) {
  cat("Classic criterion (rows) against the retired hodograph label (columns), steps:\n")
  print(table(classic = step$gaitObjective, hodograph = step$gaitHodo, useNA = "ifany"))
}

if ("outlier" %in% names(step_raw))
  cat(sprintf("Outlier filter (freq/length): removed %d/%d steps, %d/%d strides.\n",
              sum(step_raw$outlier == 1, na.rm = TRUE), nrow(step_raw),
              sum(stride_raw$outlier == 1, na.rm = TRUE), nrow(stride_raw)))

# Reconstruction quality by gait. The Methods report the root-mean-square difference
# between the force-derived and marker-based vertical CoM position per gait; that per-stride
# RMS is stored as driftRMS_vert_mm. Emit its mean (and median) over the retained strides in
# each objective gait, so the reported values are reproducible with an explicit grain rather
# than asserted in prose only. Do not write the current values into this comment. Read them from
# the emitted file.
if (all(c("gaitObjective", "driftRMS_vert_mm") %in% names(stride))) {
  if (!dir.exists("output")) dir.create("output")
  drift_by_gait <- stride %>%
    filter(!is.na(gaitObjective), is.finite(driftRMS_vert_mm)) %>%
    group_by(gait = gaitObjective) %>%
    summarise(n_strides = dplyr::n(),
              drift_mean_mm   = sprintf("%.1f", mean(driftRMS_vert_mm)),
              drift_median_mm = sprintf("%.1f", median(driftRMS_vert_mm)), .groups = "drop")
  readr::write_csv(drift_by_gait, "output/reconstruction_drift_by_gait.csv")
  cat("Reconstruction vertical-drift RMS by gait (per-stride RMS, mean over strides, mm):\n")
  print(as.data.frame(drift_by_gait))
}
if ("gaitObjective" %in% names(stride))
  print(table(hand = stride$gaitHand, objective = stride$gaitObjective))
cat(sprintf("Steadiness (strides): energy-grade (primary) %.0f%% steady vs previous energy-fraction %.0f%% vs fore-aft %.0f%% steady; grade-vs-previous agreement %.0f%%.\n",
            100*mean(stride$accClass_eg=="steady"), 100*mean(stride$accClass_en=="steady"),
            100*mean(stride$accClass_fa=="steady"), 100*mean(stride$accClass_eg==stride$accClass_en)))
