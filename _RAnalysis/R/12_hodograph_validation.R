# 12_hodograph_validation.R — the hodograph as a PROPOSED gait criterion, tested against the
# classification rather than assumed by it.
#
# Gait here is classified from the aerial phase and the classic Cavagna energy-phase criterion
# (congruity at 50%; 02_clean.R). Neither feature is computed from the CoM velocity loop, so
# every quantity below is a comparison of two independent measures.
#
# Two results:
#   (1) The counterclockwise fraction within each gait cell. The paper's one agreement
#       statistic is the unsupervised k = 2 partition against rotation sense
#       (09_gait_continuity.R).
#   (2) The CROSSOVER SPEED of each criterion, with a confidence interval: the speed at which the
#       fitted probability of a bouncing classification reaches 0.5, so half of steps fall either
#       side. It is a population quantity, not an event within a step; nothing reverses for the
#       energy-phase criterion, whose per-step cut is congruity at 50%. Both are fitted
#       with the SAME binomial GAM used for Fig 3B, P(bouncing) ~ s(u) + s(individual), over the
#       SAME sample (every classified step), and the crossing of 0.5 is read off each fit. Both
#       conditions matter: the crossing depends on which steps enter the fit, so two criteria
#       fitted to different samples are not comparable. The interval comes from a cluster
#       bootstrap over individuals, so how close the two crossovers are can be judged rather
#       than asserted.
#
# Outputs
#   output/hodo_validation_by_gait.csv     (1) per-cell counterclockwise fraction
#   output/hodo_validation_boundary.csv    (2) crossover speed of each criterion with a
#                                          cluster-bootstrap 95% confidence interval
#   output_internal/hodo_validation_boundary_fits.csv (2) the fitted curves behind those
#                                          crossover speeds; no figure draws them
suppressPackageStartupMessages({library(dplyr); library(mgcv)})
if (!dir.exists("output")) dir.create("output")

CONG <- 50
dv <- step %>% filter(is.finite(hodoArea), is.finite(congruity), is.finite(meanSpeed_n),
                      !is.na(gaitObjective)) %>%
  mutate(pend_hodo    = hodoArea > 0,        # counterclockwise = pendular sense
         pend_classic = congruity < CONG,
         # The footfall criterion, for contrast with the two mechanical ones. Taken as the
         # aerial phase rather than duty factor below 0.5: the two coincide, both are defined
         # for every step in the sample, but the aerial phase is read from the summed vertical
         # force while duty factor comes from the foot markers, which fail on a few steps
         # (a small minority below 0.25, plus those with a finite duty factor where contact
         # time is not defined).
         grounded     = hasFlight == 0)

# ---- (1) counterclockwise fraction within each gait cell -----------------------------------
n <- nrow(dv)
by_cell <- function(g, nm) dv %>% group_by(cell = as.character(.data[[g]])) %>%
  summarise(grain = nm, n_steps = dplyr::n(),
            pct_counterclockwise = round(100 * mean(pend_hodo), 1),
            .groups = "drop")
# The overall counterclockwise fraction is quoted under Statistics, as the reason the three
# agreement measures are reported together: a rare label inflates raw agreement.
all_cells <- data.frame(cell = "all steps", grain = "overall", n_steps = n,
                        pct_counterclockwise = round(100 * mean(dv$pend_hodo), 1))
readr::write_csv(bind_rows(all_cells,
                           by_cell("gaitObjective", "three-gait"),
                           by_cell("gait4", "four-cell")),
                 "output/hodo_validation_by_gait.csv")

# ---- (2) crossover speed of each criterion -------------------------------------------------
# Same model as Fig 3B. The random-effect smooth is excluded from the prediction so the curve
# is the population fit; the crossing is read on a fine grid over the observed speed range
# (lib_gam_crossing.R), and every crossing found is reported so a non-monotone fit cannot be read
# as one.
if (!exists("gam_reversal")) source("R/lib_gam_crossing.R")
cross_of <- function(dat, ybin, lab, sample_lab) {
  f <- gam_reversal(dat, ybin)
  list(row = data.frame(criterion = lab, sample = sample_lab, n_steps = f$n_obs,
                        n_crossings = f$n_crossings,
                        u_crossover = round(f$u_cross, 3),
                        Froude_crossover = round(f$froude_cross, 3),
                        p_at_slowest = round(f$p_slowest, 3),
                        p_at_fastest = round(f$p_fastest, 3)),
       fit = data.frame(criterion = lab, sample = sample_lab,
                        u = f$grid$u, p_bouncing = f$grid$p_bouncing))
}
# Both criteria are fitted over the SAME sample, every classified step, which is the sample
# Fig 3B plots. Reporting one criterion on all steps and the other on a subset is not a
# comparison, because the reversal speed depends on which steps are in the fit.
fits <- list(
  cross_of(dv, as.integer(!dv$pend_hodo),    "hodograph rotation sense", "all steps"),
  cross_of(dv, as.integer(!dv$pend_classic), "classic energy phase",     "all steps"),
  cross_of(dv, as.integer(!dv$grounded),     "aerial phase (footfall)",  "all steps"))

# ---- confidence interval on each crossover speed --------------------------------------------
# Cluster bootstrap over the individuals: steps within a bird are not independent, so the
# resampling unit is the individual, not the step. A drawn individual becomes a fresh
# random-effect level, so a bird drawn twice contributes two levels rather than one doubled one.
NBOOT <- 500
BSEED <- 7
boot_crossings <- function(dat, which_y, nboot = NBOOT, seed = BSEED) {
  set.seed(seed)
  ids <- unique(as.character(dat$subjectID))
  idx_by_id <- split(seq_len(nrow(dat)), as.character(dat$subjectID))
  out <- rep(NA_real_, nboot)
  for (b in seq_len(nboot)) {
    pick <- sample(ids, length(ids), replace = TRUE)
    idx  <- unlist(idx_by_id[pick], use.names = FALSE)
    lev  <- rep(seq_along(pick), lengths(idx_by_id[pick]))
    s <- dat[idx, ]; s$subjectID <- factor(lev)
    y <- switch(which_y,
                hodo    = as.integer(!s$pend_hodo),
                classic = as.integer(!s$pend_classic),
                aerial  = as.integer(!s$grounded))
    f <- try(gam_reversal(s, y), silent = TRUE)
    if (!inherits(f, "try-error")) out[b] <- f$u_cross
  }
  out
}
bh <- boot_crossings(dv, "hodo")
bc <- boot_crossings(dv, "classic")
ba <- boot_crossings(dv, "aerial")
ci <- function(v) stats::quantile(v, c(0.025, 0.975), na.rm = TRUE)
boundary <- bind_rows(lapply(fits, `[[`, "row")) %>%
  mutate(u_ci_lo    = round(c(ci(bh)[1], ci(bc)[1], ci(ba)[1]), 3),
         u_ci_hi    = round(c(ci(bh)[2], ci(bc)[2], ci(ba)[2]), 3),
         Froude_ci_lo = round(c(ci(bh)[1], ci(bc)[1], ci(ba)[1])^2, 3),
         Froude_ci_hi = round(c(ci(bh)[2], ci(bc)[2], ci(ba)[2])^2, 3),
         n_boot_ok  = c(sum(is.finite(bh)), sum(is.finite(bc)), sum(is.finite(ba))),
         n_boot     = NBOOT)
readr::write_csv(boundary, "output/hodo_validation_boundary.csv")
readr::write_csv(bind_rows(lapply(fits, `[[`, "fit")), "output_internal/hodo_validation_boundary_fits.csv")

cat(sprintf("12_hodograph_validation: counterclockwise fraction by gait cell over %d steps.\n", n))
print(as.data.frame(readr::read_csv("output/hodo_validation_by_gait.csv", show_col_types = FALSE)))
cat(sprintf("Crossover speed, %d-replicate cluster bootstrap over %d individuals:\n",
            NBOOT, dplyr::n_distinct(dv$subjectID)))
print(as.data.frame(boundary))
