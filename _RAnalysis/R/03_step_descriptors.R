# 03_step_descriptors.R — two derived per-step descriptors added to `step` after the QC gate.
#
# (1) Collision angle in degrees. The pipeline computes collisionAngle in radians, because the
#     force- and velocity-weighted collision angle in RADIANS is what is theoretically equal to
#     the dimensionless mechanical cost of transport (Lee et al. 2011, 2013). Degrees are easier
#     to read in a table, so both are carried: collisionAngle (rad) for the cost-of-transport
#     correspondence, collisionAngle_deg for reporting.
#
# (2) Total CoM mechanical energy at touchdown, dimensionless as E/(m g L0):
#       Ehat_TD = 1/2 uhat_TD^2 + hhat_TD,
#       uhat_TD = |v_TD| / sqrt(g L0),   hhat_TD = h_TD / L0,
#     with h_TD the CoM height above the contacting foot at touchdown, the same virtual-leg
#     convention as legLen_TD. It is computed from artifacts the production pipeline already
#     writes, so no MATLAB re-run is needed:
#       - the touchdown velocity comes from the first sample (pct = 0) of each step's cycle
#         trace in data/cycleTracesStep.csv, already normalized by sqrt(g L0). The kinetic term
#         uses the SAGITTAL CoM velocity (fore-aft and vertical), because mediolateral velocity
#         is small and every other measure in this analysis is sagittal. The cross-check below
#         records how far the sagittal form departs from the fore-aft-only one.
#       - the touchdown height comes from the virtual leg, h_TD = legLen_TD sin(alpha_TD), with
#         legLen_TD_n the touchdown leg length in units of L0 and legAngle_TD the leg angle from
#         horizontal in degrees.
#
# Plausibility gate on the touchdown height. Foot markers drop out or gap-fill to implausible
# values in a small minority of steps, and legLen_TD carries that contamination, so steps whose
# reconstructed CoM height above the foot falls outside 0.4 to 1.6 L0 are set to NA rather than
# allowed to dominate the descriptor. Ungated, the marker garbage swamps the distribution.
#
# Emits output/touchdown_energy_check.csv: the gate count, the agreement between the
# cycle-trace fore-aft touchdown velocity and the production speed_TD_n column (the recomputation
# audit for this descriptor), and the sagittal-versus-fore-aft-only difference.
suppressPackageStartupMessages({library(dplyr)})
if (!dir.exists("output")) dir.create("output")

H_TD_LIMITS <- c(0.4, 1.6)   # CoM height above the contacting foot, in units of L0

step <- step %>% mutate(collisionAngle_deg = collisionAngle * 180 / pi)

ctf <- "data/cycleTracesStep.csv"
if (!file.exists(ctf)) {
  cat("03_step_descriptors: cycleTracesStep.csv not found; touchdown energy not computed.\n")
  step$E_TD_n <- NA_real_
} else {
  td <- readr::read_csv(ctf, show_col_types = FALSE,
                        col_select = c(boutID, strideIndex, stepInStride, pct, vfa_n, vvert_n)) %>%
    group_by(boutID, strideIndex, stepInStride) %>%
    slice_min(pct, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(boutID, stepIndex = strideIndex + stepInStride - 1,
              vfa_TD_n = vfa_n, vvert_TD_n = vvert_n,
              u_TD_sag_n = sqrt(vfa_n^2 + vvert_n^2))

  step <- step %>% left_join(td, by = c("boutID", "stepIndex")) %>%
    mutate(hTD_raw_n = legLen_TD_n * sin(legAngle_TD * pi / 180),
           hTD_n     = ifelse(hTD_raw_n > H_TD_LIMITS[1] & hTD_raw_n < H_TD_LIMITS[2],
                              hTD_raw_n, NA_real_),
           E_TD_n    = 0.5 * u_TD_sag_n^2 + hTD_n)

  chk <- data.frame(
    quantity = c("steps with a touchdown cycle trace",
                 "steps gated out by the touchdown-height bound (0.4 to 1.6 L0)",
                 "steps with a finite touchdown energy",
                 "max |cycle-trace fore-aft touchdown speed - production speed_TD_n|",
                 "median (sagittal - fore-aft-only) touchdown speed, dimensionless",
                 "median touchdown height above the foot, L0"),
    value = c(sum(is.finite(step$u_TD_sag_n)),
              sum(is.finite(step$hTD_raw_n) & !is.finite(step$hTD_n)),
              sum(is.finite(step$E_TD_n)),
              round(max(abs(abs(step$vfa_TD_n) - abs(step$speed_TD_n)), na.rm = TRUE), 6),
              round(median(step$u_TD_sag_n - abs(step$vfa_TD_n), na.rm = TRUE), 4),
              round(median(step$hTD_n, na.rm = TRUE), 3)))
  readr::write_csv(chk, "output/touchdown_energy_check.csv")
  cat("03_step_descriptors: touchdown energy and collision angle in degrees added.\n")
  print(chk)
}
