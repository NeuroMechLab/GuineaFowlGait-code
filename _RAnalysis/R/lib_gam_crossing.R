# lib_gam_crossing.R — ONE definition of the crossover-speed model, shared by Fig 3
# (13_fig_rotation_sense.R) and the criterion comparison (12_hodograph_validation.R).
#
# Both scripts report the CROSSOVER SPEED, the dimensionless speed at which the fitted population
# probability of a bouncing classification reaches 0.5, so that half of steps fall either side. It
# is a property of the population of steps, not an event within a step. The
# manuscript quotes one number for it, so the model lives here rather than being written twice.
# Fitting it separately in each script let the spline basis, the grid and the speed range drift
# apart, and the two scripts then reported different values for the same quantity.
#
# Model: P(bouncing) ~ s(u, k = 6) + s(individual, bs = "re"), binomial, REML. The smoothing
# penalty is estimated from the data, so the speed term is free to shrink towards a straight line
# where that is what the data support; the returned edf_u says how far it did. The reported curve
# excludes the random-effect term, so it is the population fit. The crossing of 0.5 is
# interpolated linearly between the bracketing grid points, and EVERY crossing is counted, so a
# non-monotone fit cannot be reported as a single transition.
suppressPackageStartupMessages({library(dplyr)})

gam_reversal <- function(dat, ybin, u = "meanSpeed_n", group = "subjectID",
                         k = 6, ngrid = 400) {
  stopifnot(requireNamespace("mgcv", quietly = TRUE))
  wf <- data.frame(.y = as.integer(ybin), .u = dat[[u]], .g = factor(dat[[group]]))
  wf <- wf[is.finite(wf$.y) & is.finite(wf$.u), ]
  m <- mgcv::gam(.y ~ s(.u, k = k) + s(.g, bs = "re"), data = wf,
                 family = binomial, method = "REML")
  gu <- data.frame(.u = seq(min(wf$.u), max(wf$.u), length.out = ngrid), .g = wf$.g[1])
  gu$p <- as.numeric(predict(m, gu, type = "response", exclude = "s(.g)"))
  # Crossings of 0.5. A grid point landing exactly on 0.5 gives sign() == 0, and a raw
  # diff(sign(...)) then counts the single crossing on either side of it as two. n_crossings is
  # what says whether the fit is monotone, so carry the last non-zero sign forward before
  # differencing and count only true changes of side.
  sg <- sign(gu$p - 0.5)
  if (sg[1] == 0) sg[1] <- if (any(sg != 0)) sg[which(sg != 0)[1]] else 1
  for (i in seq_along(sg)[-1]) if (sg[i] == 0) sg[i] <- sg[i - 1]
  xs <- which(diff(sg) != 0)
  ucross <- if (length(xs))
    stats::approx(gu$p[xs[1]:(xs[1]+1)], gu$.u[xs[1]:(xs[1]+1)], 0.5)$y else NA_real_
  # Effective degrees of freedom of the speed smooth. 1 is a straight line on the logit scale,
  # so this is the positive evidence for how much curvature the fit carries.
  st <- summary(m)$s.table
  edf_u <- if ("s(.u)" %in% rownames(st)) unname(st["s(.u)", "edf"]) else NA_real_
  list(model = m, n_obs = nrow(wf), n_crossings = length(xs),
       u_cross = ucross, froude_cross = ucross^2, edf_u = edf_u, edf_u_max = k - 1,
       p_slowest = gu$p[1], p_fastest = gu$p[nrow(gu)],
       grid = data.frame(u = gu$.u, p_bouncing = gu$p))
}
