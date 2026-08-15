# 05_mean_traces.R — group-mean step-cycle traces, computed HERE in R from the per-step cycle
# traces, grouped by the gait label this pipeline assigns.
#
# The gait label lives in R (02_clean.R), so anything grouped by gait is grouped in R. Averaging
# the curves and computing the annotations printed beside them in one place is what keeps both
# describing the same set of steps, and it makes the figure follow the classifier automatically.
#
# Method. Each step contributes a 100-point trace on a percent-of-cycle axis
# (data/cycleTracesStep.csv). Within a group we average across steps at each percent and form a
# Student-t 95% confidence interval on that mean, then express the x axis in milliseconds by
# scaling percent by the group's MEAN cycle duration. Averaging on the percent axis and then
# restoring real time keeps the group's mean duration without letting a long or short cycle
# smear the landmarks, which is what averaging on an absolute time axis would do.
#
# Groupings emitted, matching what the figures consume:
#   classGait     steadiness class x gait
#   classGaitSpd  steadiness class x gait x speed half (slow/fast, split at the gait's median
#                 Froude over all classified steps, so the split does not move between classes)
#
# Output: meanTraces (in memory) and output/meanTraces_R.csv, with the same columns the figures
# expect: cycleType, grouping, grp, frbin, channel, xUnit, t_ms, mean, hi, lo, n
suppressPackageStartupMessages({library(dplyr); library(tidyr)})
if (!dir.exists("output")) dir.create("output")

CHANNELS <- c("Fz_BW", "Ffa_BW", "KE_n", "PE_n", "vfa_n", "vvert_n")
ctf <- "data/cycleTracesStep.csv"

if (!file.exists(ctf)) {
  cat("05_mean_traces: cycleTracesStep.csv not found; mean traces not built.\n")
  meanTraces <- NULL
} else {
  # per-step key -> gait, steadiness class and speed half
  spd_thresh <- step %>% filter(!is.na(gaitObjective)) %>%
    group_by(gait = as.character(gaitObjective)) %>%
    summarise(medianFroude = median(Froude, na.rm = TRUE), .groups = "drop")
  readr::write_csv(spd_thresh, "output/speedbin_thresholds_R.csv")

  keys <- step %>%
    filter(!is.na(gaitObjective), !is.na(accClass)) %>%
    transmute(boutID, stepIndex, gait = as.character(gaitObjective),
              grp = as.character(accClass), Froude) %>%
    left_join(spd_thresh, by = "gait") %>%
    mutate(spd = ifelse(Froude >= medianFroude, "fast", "slow")) %>%
    select(-Froude, -medianFroude)

  tr <- readr::read_csv(ctf, show_col_types = FALSE) %>%
    mutate(stepIndex = strideIndex + stepInStride - 1) %>%
    inner_join(keys, by = c("boutID", "stepIndex")) %>%
    select(boutID, stepIndex, gait, grp, spd, pct, cycleDur_s, all_of(CHANNELS)) %>%
    pivot_longer(all_of(CHANNELS), names_to = "channel", values_to = "value")

  # mean duration per group, used only to put the averaged trace back on a time axis
  dur_of <- function(df, keycols) df %>%
    distinct(across(all_of(c(keycols, "boutID", "stepIndex"))), cycleDur_s) %>%
    group_by(across(all_of(keycols))) %>%
    summarise(dur_s = mean(cycleDur_s, na.rm = TRUE), .groups = "drop")

  aggregate_by <- function(keycols, grouping_lab, frbin_expr) {
    dur <- dur_of(tr, keycols)
    tr %>%
      group_by(across(all_of(c(keycols, "channel", "pct")))) %>%
      summarise(mean = mean(value, na.rm = TRUE),
                sd   = stats::sd(value, na.rm = TRUE),
                n    = sum(is.finite(value)), .groups = "drop") %>%
      # Student-t 95% interval on the group mean at each percent of the cycle
      mutate(hw = ifelse(n > 1, stats::qt(0.975, n - 1) * sd / sqrt(n), NA_real_),
             hi = mean + hw, lo = mean - hw) %>%
      left_join(dur, by = keycols) %>%
      mutate(cycleType = "step", grouping = grouping_lab, xUnit = "ms",
             t_ms = pct / 100 * dur_s * 1000,
             frbin = frbin_expr(pick(everything()))) %>%
      select(cycleType, grouping, grp, frbin, channel, xUnit, t_ms, mean, hi, lo, n)
  }

  meanTraces <- bind_rows(
    aggregate_by(c("grp", "gait"),        "classGait",    function(d) d$gait),
    aggregate_by(c("grp", "gait", "spd"), "classGaitSpd", function(d) paste(d$gait, d$spd, sep = "|"))
  )
  readr::write_csv(meanTraces, "output/meanTraces_R.csv")

  # Per-cell annotation over the EXACT steps behind each mean trace, so the strip text and the
  # curve describe the same set. Emitted for the figure scripts and for the captions.
  ann_of <- function(keycols, frbin_fun) {
    st <- step %>% filter(!is.na(gaitObjective), !is.na(accClass)) %>%
      mutate(gait = as.character(gaitObjective), grp = as.character(accClass)) %>%
      left_join(spd_thresh, by = "gait") %>%
      mutate(spd = ifelse(Froude >= medianFroude, "fast", "slow"))
    st %>% group_by(across(all_of(keycols))) %>%
      summarise(n = dplyr::n(),
                Fr_mean = mean(Froude), Fr_min = min(Froude), Fr_max = max(Froude),
                u_mean = mean(meanSpeed_n), u_min = min(meanSpeed_n), u_max = max(meanSpeed_n),
                v_mean = mean(meanSpeed), v_min = min(meanSpeed), v_max = max(meanSpeed),
                .groups = "drop") %>%
      mutate(frbin = frbin_fun(pick(everything())))
  }
  meanTracesAnn <- bind_rows(
    ann_of(c("grp", "gait"),        function(d) d$gait) %>% mutate(grouping = "classGait"),
    ann_of(c("grp", "gait", "spd"), function(d) paste(d$gait, d$spd, sep = "|")) %>%
      mutate(grouping = "classGaitSpd"))
  readr::write_csv(meanTracesAnn, "output/meanTraces_annot_R.csv")

  # Timing of the peak of the mean vertical force within contact, as a fraction of contact
  # duration. The Discussion attributes the early-peaking avian force profile to leg damping, so
  # the observation behind that has to be a measured number rather than a look at the figure.
  peak_timing <- function(grouping_lab) {
    meanTraces %>% filter(grouping == grouping_lab, grp == "steady", channel == "Fz_BW") %>%
      group_by(grouping, frbin) %>% arrange(t_ms, .by_group = TRUE) %>%
      summarise(peak_BW = round(max(mean), 2),
                pct_of_contact = { thr <- 0.05 * max(mean); on <- which(mean > thr)
                  round(100 * (t_ms[which.max(mean)] - t_ms[min(on)]) /
                        (t_ms[max(on)] - t_ms[min(on)])) },
                .groups = "drop")
  }
  fz_peak <- bind_rows(peak_timing("classGait"), peak_timing("classGaitSpd"))
  readr::write_csv(fz_peak, "output/vertical_force_peak_timing.csv")
  cat("05_mean_traces: peak of the mean vertical force, percent of contact:\n")
  print(as.data.frame(fz_peak))

  cat(sprintf(paste0("05_mean_traces: %d steps entered the group means; %d groups over %d ",
                     "channels; classGaitSpd cells: %s\n"),
              dplyr::n_distinct(paste(tr$boutID, tr$stepIndex)),
              dplyr::n_distinct(paste(meanTraces$grouping, meanTraces$grp, meanTraces$frbin)),
              length(CHANNELS),
              paste(sort(unique(meanTraces$frbin[meanTraces$grouping == "classGaitSpd"])),
                    collapse = ", ")))
  print(as.data.frame(meanTracesAnn %>% filter(grouping == "classGaitSpd", grp == "steady") %>%
                        select(gait, spd, n, Fr_mean, u_mean)))
}
