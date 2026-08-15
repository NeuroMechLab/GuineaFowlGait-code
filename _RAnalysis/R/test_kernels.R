#!/usr/bin/env Rscript
# test_kernels.R - known-answer tests for the computations a reader cannot check by eye.
#
# Not sourced by run_all.R and it writes nothing: run it by hand after changing any of the
# kernels below. Every case has an answer derived independently of the implementation, either
# analytically or by construction, and each edge case is one the pipeline has actually met.
#
#   Rscript R/test_kernels.R        (from _RAnalysis)
#
# Covers: the shoelace signed area and its sign convention, the Blum 2009 method-C leg-compression
# inversion and its guards, the merged and split mis-cut flags, and the crossover reader on a
# monotone curve and on a non-monotone one that must report two crossings.
suppressPackageStartupMessages({library(dplyr)})
ok <- 0L; bad <- 0L
chk <- function(label, got, want, tol = 1e-9) {
  pass <- isTRUE(all.equal(got, want, tolerance = tol))
  if (pass) ok <<- ok + 1L else bad <<- bad + 1L
  cat(sprintf("%-4s %-62s got %s want %s\n", if (pass) "PASS" else "FAIL", label,
              paste(format(got, digits = 6), collapse = ","),
              paste(format(want, digits = 6), collapse = ",")))
}

# ---- 1. shoelace signed area -----------------------------------------------------------------
# A unit circle traversed counterclockwise encloses +pi; reversing the traversal negates it. The
# pipeline's convention is positive = counterclockwise = pendular, so the sign is the result.
shoelace <- function(x, y) { n <- length(x); j <- c(2:n, 1)
  0.5 * sum(x * y[j] - x[j] * y) }
th <- seq(0, 2 * pi, length.out = 2001)[-2001]
chk("shoelace, unit circle counterclockwise", shoelace(cos(th), sin(th)), pi, tol = 1e-5)
chk("shoelace, unit circle clockwise",        shoelace(cos(-th), sin(-th)), -pi, tol = 1e-5)
# a 3-4-5 right triangle, area 6 exactly
chk("shoelace, 3-4-5 triangle", shoelace(c(0, 3, 0), c(0, 0, 4)), 6)
# translation invariance: the signed area must not move with the origin, which is why the
# measure survives a constant error in the vertical-velocity integration constant
chk("shoelace, translation invariant", shoelace(cos(th) + 17, sin(th) - 4), pi, tol = 1e-5)
# two equal, opposite lobes net to zero, which is why a self-intersecting loop is reported by the
# sense that dominates BY AREA rather than being given a hybrid category
chk("shoelace, equal opposite lobes net to zero",
    shoelace(c(cos(th), rev(cos(th))), c(sin(th), rev(sin(th)))), 0, tol = 1e-5)
chk("shoelace, unequal opposite lobes keep the larger sense",
    sign(shoelace(c(2 * cos(th), rev(cos(th))), c(2 * sin(th), rev(sin(th))))), 1)

# ---- 2. method-C leg compression and stiffness ------------------------------------------------
# dL = L0 + (Fmax/m)(tc/pi)^2 - (g/8)tc^2 - L0 sin(alpha); khat = (Fmax/(m g)) / (dL/L0).
g <- 9.81
methodC <- function(L0, m, Fmax, tc, alpha_deg) {
  dL <- L0 + (Fmax / m) * (tc / pi)^2 - (g / 8) * tc^2 - L0 * sin(alpha_deg * pi / 180)
  frac <- dL / L0
  khat <- if (is.finite(frac) && frac > 0 && frac < 0.40) (Fmax / (m * g)) / frac else NA_real_
  c(frac = frac, khat = khat)
}
# hand-computed case: alpha = 90 deg makes sin(alpha) = 1, so the L0 terms cancel exactly and
# dL = (Fmax/m)(tc/pi)^2 - (g/8)tc^2, independent of L0.
m <- 1.4; Fmax <- 2 * m * g; tc <- 0.20; L0 <- 0.21
dL_expect <- (Fmax / m) * (tc / pi)^2 - (g / 8) * tc^2
r <- methodC(L0, m, Fmax, tc, 90)
chk("method C, alpha = 90 deg, L0 terms cancel", unname(r["frac"]), dL_expect / L0)
chk("method C, khat = (Fmax/BW)/(dL/L0)", unname(r["khat"]), 2 / (dL_expect / L0))
# the angle enters only through sin, so alpha and 180 - alpha must agree. This is why a touchdown
# angle of 122.8 deg and its acute complement give the same stiffness.
chk("method C, sin symmetry about 90 deg",
    unname(methodC(L0, m, Fmax, tc, 122.8)["khat"]),
    unname(methodC(L0, m, Fmax, tc, 57.2)["khat"]))
# guards: a leg that extends over stance gives dL <= 0 and must return NA, not a negative
# stiffness. Vaulting steps do this, which is why walking loses the most steps.
chk("method C, non-physical dL <= 0 returns NA",
    unname(methodC(L0, m, 0.2 * m * g, 0.05, 20)["khat"]), NA_real_)
chk("method C, compression above 0.4 L0 returns NA",
    unname(methodC(L0, m, 12 * m * g, 0.35, 20)["khat"]), NA_real_)

# ---- 3. merged and split mis-cut flags --------------------------------------------------------
# Forward speed = step length x step frequency, so a merged cycle halves the frequency and a split
# cycle doubles it while the speed is unchanged. The flag is |log(f / f_local_median)| > |log(0.65)|.
DBL_FRAC <- 0.65
flag_double <- function(freq, med) abs(log(freq / med)) > abs(log(DBL_FRAC))
med <- 4.0
chk("mis-cut, a clean cycle is not flagged",       flag_double(4.10, med), FALSE)
chk("mis-cut, a merged (half-frequency) cycle",    flag_double(2.00, med), TRUE)
chk("mis-cut, a split (double-frequency) cycle",   flag_double(8.00, med), TRUE)
# symmetry: the criterion must catch both directions equally, which a one-sided test would not
chk("mis-cut, symmetric in log frequency",
    flag_double(med * 0.64, med), flag_double(med / 0.64, med))

# ---- 4. crossover reader -----------------------------------------------------------------------
# Read every crossing of 0.5, so a non-monotone fit cannot be reported as a single transition.
# same logic as gam_reversal in lib_gam_crossing.R: carry the last non-zero sign forward so a
# grid point landing exactly on 0.5 is not counted as two crossings
crossings <- function(u, p, lev = 0.5) {
  sg <- sign(p - lev)
  if (sg[1] == 0) sg[1] <- if (any(sg != 0)) sg[which(sg != 0)[1]] else 1
  for (i in seq_along(sg)[-1]) if (sg[i] == 0) sg[i] <- sg[i - 1]
  i <- which(diff(sg) != 0)
  if (!length(i)) return(numeric(0))
  u[i] + (lev - p[i]) * (u[i + 1] - u[i]) / (p[i + 1] - p[i])
}
u <- seq(0, 2, by = 0.001)
chk("crossover, logistic crossing at u = 0.8",
    crossings(u, 1 / (1 + exp(-12 * (u - 0.8)))), 0.8, tol = 1e-3)
# A curve that starts below 0.5, rises through it and falls back must report TWO crossings, at
# u = 0.25 and u = 1.25 by construction. Reporting one would hide exactly the kind of shape the
# leg-stiffness result turned on. (A curve that merely TOUCHES 0.5 at its ends does not cross
# there, which is why the phase offset is needed to make this a genuine double crossing.)
p2 <- 0.5 + 0.4 * sin(pi * (u - 0.25))
chk("crossover, non-monotone curve reports two", length(crossings(u, p2)), 2L)
chk("crossover, non-monotone crossings located", round(crossings(u, p2), 2), c(0.25, 1.25))
chk("crossover, a curve that never reaches 0.5 reports none",
    length(crossings(u, rep(0.2, length(u)))), 0L)
# a grid point landing exactly on 0.5 must still be ONE crossing, not two
chk("crossover, an exact 0.5 grid hit counts once",
    length(crossings(c(0, 1, 2), c(0.2, 0.5, 0.8))), 1L)

cat(sprintf("\n%d passed, %d failed\n", ok, bad))
if (bad > 0) quit(status = 1)
