# lib_format.R — shared summary-cell formatters for the paper tables (19_table_gait_summary.R,
# 20_table_gait4_summary.R), so every table in the manuscript uses one format.
#
# Every distributional cell is median (interquartile range). Several of the reported measures are
# skewed or bounded (leg stiffness has a right tail; duty factor is bounded near 0.5;
# recovery is bounded below), so a median and IQR describe a group without a mean being pulled by
# a tail. No table reports min-max bounds.
#
# Rounding is set by the data rather than chosen per measure. For a measure the decimal count is
# the number of places needed to show the NARROWEST interquartile width among the compared groups
# to two significant figures,
#   decimals = max(0, 1 - floor(log10(w_min))).
# Reporting a median more finely than the spread supports would imply a precision the data do not
# have. The narrowest rather than the typical width sets the count because a measure whose spread
# differs between groups by an order of magnitude (Froude number, whose walking interquartile
# width is a tenth of the aerial-running one) would otherwise be rounded until the tightest group
# collapsed to a single repeated digit. Erring toward one extra place in the widest group is the
# safer of the two errors.

iqr_digits <- function(x_by_group) {
  w <- vapply(x_by_group, function(x) { x <- x[is.finite(x)]
    if (length(x) < 4) return(NA_real_)
    unname(diff(stats::quantile(x, c(0.25, 0.75)))) }, numeric(1))
  w <- w[is.finite(w) & w > 0]
  if (!length(w)) return(2L)
  as.integer(max(0, min(5, 1 - floor(log10(min(w))))))
}

medIQR <- function(x, d) { x <- x[is.finite(x)]
  if (!length(x)) return("-")
  q <- stats::quantile(x, c(0.25, 0.5, 0.75))
  sprintf("%.*f (%.*f to %.*f)", d, q[2], d, q[1], d, q[3]) }

# Percentage of a signed quantity above and below zero, for the rotation-sense consistency row.
rot <- function(a) { a <- a[is.finite(a)]
  if (!length(a)) return("-")
  sprintf("%.0f / %.0f", 100 * mean(a > 0), 100 * mean(a < 0)) }
