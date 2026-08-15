# lib_phase_register.R — Hilbert continuous-phase registration for step/stride cycles.
#
# Percent-of-cycle averaging aligns only the cycle endpoints, so variation in the timing of the
# force peak within a step, and in how a stride divides between its two steps, smears and lowers
# the averaged peaks.
# Instead this derives ONE continuous phase for the whole stride from the analytic
# signal (Hilbert transform) of the vertical CoM velocity, a once-per-step
# oscillator for every gait, and resamples the channels at equal phase. Phase
# advances continuously across the step transition (no step-join discontinuity),
# and peaks/troughs align because they recur at a consistent phase of the
# vertical-velocity oscillation. The step transition lands at a data-determined
# phase (returned as `bnd`), not a forced landmark. Used by 14_fig_hodographs.R
# (steady strides) and 17_fig_transitions.R (transition strides).
suppressPackageStartupMessages({library(dplyr)})

PR_CHAN <- c("Fz_BW","Ffa_BW","KE_n","PE_n","vfa_n","vvert_n")

.pr_analytic <- function(x) {                    # analytic signal via FFT (discrete Hilbert)
  n <- length(x); X <- fft(x); h <- numeric(n)
  if (n %% 2 == 0) { h[c(1, n/2 + 1)] <- 1; h[2:(n/2)] <- 2 }
  else            { h[1] <- 1; h[2:((n + 1)/2)] <- 2 }
  fft(X * h, inverse = TRUE) / n
}
.pr_unwrap <- function(p) {                      # phase unwrap (no base-R builtin)
  d <- diff(p); d <- d - 2*pi*round(d/(2*pi)); cumsum(c(p[1], d))
}
# continuous stride phase (0..1) from a once-per-step oscillator; the signal is
# tiled x3 (stride endpoints are both touchdown, so it is ~periodic) and only the
# middle third kept, to suppress FFT edge effects. The instantaneous frequency
# (phase increment) is lightly smoothed and floored at a small positive value, so
# phase rises strictly and monotonically: this removes the phase-slip artifact
# (where low oscillation amplitude stalls the raw phase, then a compensating jump
# distorts the resampling) without a hard cummax plateau.
.pr_stride_phase <- function(x) {
  n <- length(x); xc <- x - mean(x)
  ph <- .pr_unwrap(Arg(.pr_analytic(rep(xc, 3))))[(n + 1):(2*n)]
  d <- diff(ph)
  k <- max(3, round(n / 40))
  ds <- as.numeric(stats::filter(d, rep(1/k, k), sides = 2))
  ds[is.na(ds)] <- d[is.na(ds)]                 # keep endpoints where the window overruns
  ds <- pmax(ds, 1e-6)                          # strictly increasing, no plateau
  ph <- cumsum(c(0, ds))
  (ph - ph[1]) / (ph[n] - ph[1])
}
# Hilbert-registered transition strides: one continuous phase 0..200 per stride
# (step transition at the data-determined phase `bnd`), channels resampled at
# equal phase, then averaged across strides. Returns one row per (strideID, phase)
# with the registered PR_CHAN channels, `bnd`, and `vfa_c` (fore-aft velocity
# fluctuation about the stride mean, for hodographs) + transClass. `ctStep` =
# cycleTracesStep; `trans` (or the steady-stride table) supplies strideID/transClass.
pr_registered_transitions_hilbert <- function(ctStep, trans, npts = 200) {
  key <- trans %>% mutate(strideID = paste(boutID, strideIndex)) %>% select(strideID, transClass)
  ct <- ctStep %>% mutate(strideID = paste(boutID, strideIndex)) %>% inner_join(key, by = "strideID")
  grid <- seq(0, 1, length.out = npts)
  regs <- ct %>% group_by(strideID, transClass) %>% group_modify(~{
    d <- .x %>% mutate(u = ifelse(stepInStride == 1, pct/100, 1 + pct/100)) %>%
      arrange(u) %>% distinct(u, .keep_all = TRUE)
    phn <- .pr_stride_phase(d$vvert_n)
    res <- data.frame(phase = grid * 200, bnd = stats::approx(d$u, phn, 1, rule = 2)$y * 200)
    for (ch in PR_CHAN) res[[ch]] <- stats::approx(phn, d[[ch]], grid, rule = 2)$y
    res
  }) %>% ungroup()
  regs %>% group_by(strideID) %>% mutate(vfa_c = vfa_n - mean(vfa_n)) %>% ungroup()
}
