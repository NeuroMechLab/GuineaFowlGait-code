# 06_fig_outlier_filter.R — outlier-filter validation graph (supplement). Shows the step and stride
# frequency/length distributions and the frequency-versus-speed relationship, with the flagged
# (removed) observations coloured, so a human can confirm the filter behaves sensibly before
# trusting the downstream analysis. The filter is SPEED-CONDITIONAL (a log-MAD frequency/length
# filter judged within dimensionless-speed neighbourhoods, plus a speed-conditional doubling
# filter), so there is no single global keep-box: the flagged points are the cycles that are
# anomalous FOR THEIR SPEED (the merged/split mis-cuts that fall off the frequency-speed trend).
# Reads the RAW tables (step_raw / stride_raw), which carry the authoritative 'outlier' flag.
suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(patchwork)})
if (!exists("theme_daley")) source("R/lib_theme.R")
if (!"outlier" %in% names(step_raw)) {
  cat("06_fig_outlier_filter: no outlier flag found; skipping.\n")
} else {
  ocol <- c(`0` = "grey55", `1` = "#C0398B")
  olab <- c(`0` = "kept", `1` = "removed")
  flagf <- function(df) df %>% mutate(o = factor(ifelse(outlier == 1, 1, 0), levels = c(0,1)))

  hist_panel <- function(df, var, xlab) {
    d <- flagf(df %>% filter(is.finite(.data[[var]]), .data[[var]] > 0))
    ggplot(d, aes(.data[[var]], fill = o)) +
      geom_histogram(bins = 45, colour = NA) + scale_x_log10() +
      scale_fill_manual(values = ocol, labels = olab, name = NULL, drop = FALSE) +
      labs(x = xlab, y = "count") + theme_bio()
  }
  pA <- hist_panel(step_raw,   "stepFreq",   "step frequency (Hz)")
  pB <- hist_panel(step_raw,   "stepLength", "step length (m)")
  pC <- hist_panel(stride_raw, "strideFreq", "stride frequency (Hz)")
  pD <- hist_panel(stride_raw, "strideLength", "stride length (m)")

  # frequency versus dimensionless speed: kept cycles trace a tight band; the flagged cycles
  # are the merged (~half) and split (~double) mis-cuts that depart from the speed-local band.
  speed_panel <- function(df, fvar, xlab) {
    if (!all(c("meanSpeed_n", fvar) %in% names(df))) return(NULL)
    d <- flagf(df %>% filter(is.finite(meanSpeed_n), meanSpeed_n > 0,
                             is.finite(.data[[fvar]]), .data[[fvar]] > 0))
    ggplot(d, aes(meanSpeed_n, .data[[fvar]], colour = o)) +
      geom_point(alpha = 0.55, size = 0.7) + scale_y_log10() +
      scale_colour_manual(values = ocol, labels = olab, name = NULL, drop = FALSE) +
      labs(x = "dimensionless speed u", y = xlab) + theme_bio() +
      # One kept/removed key serves the whole figure; the histogram fill carries it.
      guides(colour = "none")
  }
  pE <- speed_panel(step_raw,   "stepFreq",   "step frequency (Hz)")
  pF <- speed_panel(stride_raw, "strideFreq", "stride frequency (Hz)")

  nStepOut <- sum(step_raw$outlier == 1, na.rm = TRUE)
  nStrOut  <- sum(stride_raw$outlier == 1, na.rm = TRUE)
  # No figure title or subtitle: the number and the description live in the caption only, so
  # they cannot drift out of step with the manuscript when figures are renumbered.
  # Each panel in a column carries a different measure on its own range, so every axis is
  # labelled.
  fig <- ((pA | pB) / (pC | pD) / (pE | pF)) + plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "A")
  bio_save("SFig_OutlierFilter", fig, width = BIO_W2, height = 7.4)
  cat(sprintf("06_fig_outlier_filter: SFig_OutlierFilter (removed %d steps, %d strides; speed-conditional).\n", nStepOut, nStrOut))
}
