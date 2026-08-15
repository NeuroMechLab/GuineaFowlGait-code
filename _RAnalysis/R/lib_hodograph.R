# lib_hodograph.R — geometry of the CoM velocity loop.
#
# Every measure here reads the SAME closed polygon that the classifier's signed area reads: the
# samples of one cycle in the (fore-aft velocity fluctuation, vertical velocity) plane, closed
# from the last sample back to the first. Keeping them in one file is what stops a loop
# descriptor from being defined on a curve the signed area was not computed on, which is the way
# an annotation ends up contradicting the gait label it annotates.
#
# Callers pass x and y already normalized by sqrt(g L0); centring is done here where a measure
# needs it, because the shoelace sum is translation invariant while the turning angle is not.

# Shoelace signed area of the closed polygon. Positive is counterclockwise (pendular),
# negative clockwise (bouncing). For a self-intersecting loop this is the algebraic sum of the
# lobe areas, so the sign is the sense that dominates by area.
signed_area <- function(x, y) {
  xn <- c(x[-1], x[1]); yn <- c(y[-1], y[1])
  0.5 * sum(x * yn - xn * y)
}

# How the loop's area divides between the two senses. Each segment of the closed trace
# contributes 0.5 * (x_i y_{i+1} - x_{i+1} y_i) to the shoelace sum, the area swept by the radius
# vector, and those contributions add to the signed area itself. Splitting them by sign therefore
# decomposes the very quantity that classifies the cycle: the share above 50% IS the net sense, so
# the split cannot contradict the gait label. The share also says how decisively the loop turns one
# way, which a bare sign does not: a cycle near the crossover splits its area almost evenly while a
# fast aerial run is nearly all of one sense.
# One decimal, not zero: a loop that is 99.6% one sense must not print as 100%.
area_split <- function(x, y) {
  xc <- x - mean(x); yc <- y - mean(y)
  dA <- 0.5 * (xc * c(yc[-1], yc[1]) - c(xc[-1], xc[1]) * yc)   # sum(dA) is the signed area
  gross <- sum(abs(dA))
  c(pct_area_cw  = 100 * sum(abs(dA[dA < 0])) / gross,
    pct_area_ccw = 100 * sum(dA[dA > 0]) / gross)
}

# Share of the loop's gross area enclosed in its NET sense, which is the larger of the two shares
# above. 50% is an even split, 100% a loop that turns one way throughout.
dominant_area_share <- function(x, y) max(area_split(x, y))

# Fraction of the total turning angle of the radius vector that runs AGAINST the net sense. A
# local reversal is graded by this rather than made a third category: it matters only if it
# dominates, and if it does the signed area has already changed sign.
loop_retro_frac <- function(x, y) {
  xc <- x - mean(x); yc <- y - mean(y)
  th <- atan2(yc, xc)
  d  <- diff(c(th, th[1]))
  d  <- ((d + pi) %% (2 * pi)) - pi          # wrap each increment to (-pi, pi]
  tot <- sum(d); den <- sum(abs(d))
  retro <- if (tot >= 0) sum(abs(d[d < 0])) else sum(abs(d[d > 0]))
  if (den > 0) retro / den else NA_real_
}

# Number of times the trajectory crosses itself: a figure-of-eight returns 1, a simple loop 0.
# Two segments count as a crossing only when each strictly straddles the other's line, so shared
# endpoints do not, and consecutive segments are excluded from the comparison for the same reason.
# A crossing splits the curve into lobes of opposite sense, whose areas the signed area then sums,
# so this is a property of the drawn curve and not of the sign of its area.
#
# closed = FALSE, the default, counts crossings of the TRAVERSED path, the n - 1 segments between
# consecutive samples. A step cycle does not return to its starting velocity, so the closing
# segment from the last sample back to the first is a straight chord the bird never traversed:
# a crossing against that chord is a property of how the shoelace formula closes the polygon, not
# a crossover of the animal's velocity trajectory, so it must not be reported as one.
# closed = TRUE counts them, for the polygon the signed area is computed on.
loop_self_crossings <- function(x, y, closed = FALSE) {
  n <- length(x)
  if (n < 4) return(NA_integer_)
  if (closed) { xa <- x; ya <- y; xb <- c(x[-1], x[1]); yb <- c(y[-1], y[1]) }
  else        { xa <- x[-n]; ya <- y[-n]; xb <- x[-1]; yb <- y[-1] }
  m <- length(xa)
  cnt <- 0L
  for (i in seq_len(m - 2)) {
    # segments i and i+1 always share a vertex; on the closed polygon so do segments m and 1
    hi <- if (closed && i == 1L) m - 1L else m
    if (i + 2L > hi) next
    j  <- (i + 2L):hi
    ex <- xb[i] - xa[i]; ey <- yb[i] - ya[i]
    fx <- xb[j] - xa[j]; fy <- yb[j] - ya[j]
    d1 <- fx * (ya[i] - ya[j]) - fy * (xa[i] - xa[j])
    d2 <- fx * (yb[i] - ya[j]) - fy * (xb[i] - xa[j])
    d3 <- ex * (ya[j] - ya[i]) - ey * (xa[j] - xa[i])
    d4 <- ex * (yb[j] - ya[i]) - ey * (xb[j] - xa[i])
    cnt <- cnt + sum(d1 * d2 < 0 & d3 * d4 < 0, na.rm = TRUE)
  }
  as.integer(cnt)
}
