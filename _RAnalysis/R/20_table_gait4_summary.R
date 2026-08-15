# 20_table_gait4_summary.R — four-cell CoM-dynamics summary tables.
#
# One main table over steady steps and one supplement table over unsteady steps
# (accelerating above, decelerating below), across the four cells of the classifier's two
# qualitative features, aerial phase x pendular-vs-bouncing energy phase (gait4, 02_clean.R):
#   Walk (grounded)  grounded + pendular
#   Grounded run     grounded + bouncing
#   Aerial run       aerial   + bouncing
#   Pendular run     aerial   + pendular   (split out of the gaitObjective aerialRun class)
#
# These support Results 1-2: they give the per-cell descriptive statistics behind the
# continuum claim, and they quantify the aerial pendular cell instead of asserting it empty.
#
# The rotation-sense row is an OBSERVED result here, not a restatement of the cell
# definitions. The cells are defined by the aerial phase and by KE-PE congruity at 50%,
# neither of which uses the CoM velocity loop, so the counterclockwise fraction of each cell
# is a measured correspondence between independent quantities.
#
# PC1 and PC2 are the CoM-dynamics principal components, refitted here with the same
# features, subset and sign convention as 18_sfig_pc_steadiness.R and Table S1
# (09_gait_continuity.R) so all three quote the same component.
#
# Writes output/GaitSummary_FourCell.xlsx (tables, per-cell n, composition) plus one CSV per
# table and a markdown rendering for drafting.
suppressPackageStartupMessages({library(dplyr); library(openxlsx)})
if (!exists("theme_daley")) source("R/lib_theme.R")
if (!dir.exists("output")) dir.create("output")

# ---- formatters -------------------------------------------------------------
# median (IQR) cells with data-driven rounding, shared with 19_table_gait_summary.R so every table in the
# manuscript uses one format. Digits here are fixed over ALL steadiness classes, so the steady and
# unsteady tables share one format per measure.
if (!exists("medIQR")) source("R/lib_format.R")

# ---- PC scores, matching the PC-space supplement ----------------------------
# Must stay identical to `feats` in 09_gait_continuity.R and 18_sfig_pc_steadiness.R: the
# three places report the same components, and the hodograph signed area is not an input
# because it is the measure tested against this space in 12_hodograph_validation.R.
PCA_FEATS <- SAMPLE_DESCRIPTORS   # the gated descriptor set (04_analysis_sample.R)
dp <- step %>% filter(as.character(gaitObjective) %in% GAIT_LEVELS,
                      if_all(all_of(PCA_FEATS), is.finite), !is.na(accClass))
pcFit <- prcomp(scale(as.matrix(dp[, PCA_FEATS])))
s1 <- if (mean(pcFit$x[dp$gaitObjective == "walk", 1]) < 0) -1 else 1
s2 <- if (mean(pcFit$x[dp$accClass == "steady", 2]) >
           mean(pcFit$x[dp$accClass != "steady", 2])) -1 else 1
dp$PC1 <- pcFit$x[, 1] * s1; dp$PC2 <- pcFit$x[, 2] * s2
stepT <- step %>% left_join(dp %>% select(boutID, stepIndex, PC1, PC2),
                            by = c("boutID", "stepIndex"))

# ---- measure registry -------------------------------------------------------
# The first row gives the steady step count alongside the count over ALL steadiness classes, so
# the reader can see at once what fraction of each cell the summary rests on; the unsteady
# remainder is summarized in the supplement table. TOT_BY_CELL is the all-class denominator.
TOT_BY_CELL <- table(factor(as.character(stepT$gait4), levels = GAIT4_LEVELS))
MEASURES <- list(
  list(label = "Number of steps (this class of all)", col = NA,   kind = "count"),
  list(label = "Dimensionless speed, u",                  col = "meanSpeed_n",        kind = "medIQR"),
  list(label = "Froude number, u^2",                      col = "Froude",             kind = "medIQR"),
  list(label = "PC1 score",                               col = "PC1",                kind = "medIQR"),
  list(label = "PC2 score",                               col = "PC2",                kind = "medIQR"),
  list(label = "Hodograph rotation sense (% CCW / % CW)",  col = "hodoArea",           kind = "rot"),
  list(label = "Hodograph signed area, A/(g L0)",          col = "hodoArea_n",         kind = "medIQR"),
  list(label = "Pendular recovery (%)",                   col = "recovery",           kind = "medIQR"),
  list(label = "KE-PE congruity (%)",                     col = "congruity",          kind = "medIQR"),
  list(label = "Mechanical cost of transport",            col = "CoTmech",            kind = "medIQR"),
  list(label = "Collision angle (deg)",                   col = "collisionAngle_deg", kind = "medIQR"),
  list(label = "Touchdown mechanical energy, E/(m g L0)", col = "E_TD_n",             kind = "medIQR"),
  list(label = "Leg stiffness, k L0/(m g)",               col = "kLeg_n",             kind = "medIQR"),
  # Leg stiffness is the one measure outside the analysis gate, so its cells rest on fewer steps
  # than the column count in row 1. That n is given here rather than left to the cellN sheet,
  # because a reader takes the row-1 count to hold for every row below it.
  list(label = "Steps with a leg-stiffness fit",          col = "kLeg_n",             kind = "nfin"),
  list(label = "Duty factor",                            col = "dutyFactor",         kind = "medIQR"))

# Decimal places per measure, from the narrowest interquartile width over ALL steadiness classes
# pooled (see iqr_digits), so Table 1 and Table S4 share one format for every measure.
#
# Setting the width from the steady cells alone was tried and reverted: it takes Froude to 3 dp,
# which is the stricter reading of the convention for Table 1's own narrowest Froude width, but it
# also takes KE-PE congruity to whole percent, because the narrowest congruity width among the
# steady cells is wide. Whole-percent congruity contradicts the values Results quotes. Pooling keeps both columns readable; the cost is that Froude, a derived convenience
# column beside the dimensionless speed the paper actually reports, shows its narrowest width to
# one significant figure.
DIGITS <- vapply(MEASURES, function(m) {
  if (m$kind != "medIQR") return(0L)
  iqr_digits(lapply(GAIT4_LEVELS, function(g) stepT[[m$col]][which(stepT$gait4 == g)]))
}, integer(1))
fmt_cell <- function(m, d, dat, g) switch(m$kind,
  count  = sprintf("%d of %d", nrow(dat), TOT_BY_CELL[[g]]),
  rot    = rot(dat[[m$col]]),
  nfin   = sprintf("%d", sum(is.finite(dat[[m$col]]))),
  medIQR = medIQR(dat[[m$col]], d))

build <- function(dat) {
  cells <- lapply(GAIT4_LEVELS, function(g) dat %>% filter(gait4 == g))
  names(cells) <- GAIT4_LEVELS
  tab <- data.frame(Measure = vapply(MEASURES, function(m) m$label, ""), stringsAsFactors = FALSE)
  ntab <- tab
  for (g in GAIT4_LEVELS) {
    nm <- unname(GAIT4_LABELS[g])
    tab[[nm]]  <- vapply(seq_along(MEASURES),
                         function(i) fmt_cell(MEASURES[[i]], DIGITS[i], cells[[g]], g), "")
    ntab[[nm]] <- vapply(MEASURES, function(m) if (is.na(m$col)) nrow(cells[[g]])
                         else sum(is.finite(cells[[g]][[m$col]])), numeric(1))
  }
  list(tab = tab, n = ntab)
}

main <- build(stepT %>% filter(accClass == "steady"))
acc  <- build(stepT %>% filter(accClass == "accelerating"))
dec  <- build(stepT %>% filter(accClass == "decelerating"))
sup   <- bind_rows(data.frame(Section = "Accelerating steps", acc$tab, check.names = FALSE),
                   data.frame(Section = "Decelerating steps", dec$tab, check.names = FALSE))
sup_n <- bind_rows(data.frame(Section = "Accelerating steps", acc$n,   check.names = FALSE),
                   data.frame(Section = "Decelerating steps", dec$n,   check.names = FALSE))

# Composition of each cell by steadiness class. This is what makes the aerial pendular cell
# interpretable: not merely rare, but overwhelmingly decelerating.
comp <- stepT %>% filter(!is.na(gait4), !is.na(accClass)) %>%
  count(gait4, accClass) %>%
  tidyr::pivot_wider(names_from = accClass, values_from = n, values_fill = 0) %>%
  mutate(total = decelerating + steady + accelerating,
         pct_steady = round(100 * steady / total),
         pct_decelerating = round(100 * decelerating / total),
         gait = unname(GAIT4_LABELS[as.character(gait4)])) %>%
  select(gait, total, decelerating, steady, accelerating, pct_steady, pct_decelerating)

# Which observations each figure draws on. The gait classification uses every QC-passed step, so
# steadiness is a label rather than a filter, and the figures split accordingly: the two
# classification-evidence figures use all classified observations, and
# the descriptor and speed-relation figures use the steady subset. Emitted so the Table 1
# footnote quotes counts from a file rather than from prose.
# the Tables S1 to S3 sample: the descriptor axes those analyses require
PCA_FEATS_FIG <- PCA_FEATS
fig_subsets <- data.frame(
  figure = c("Tables S1-S3","Fig 2","Fig 3","Fig 4","Fig 5","Fig 6","Fig 7"),
  observations = c("steps","steps","steps","steps","strides","steps","steps"),
  subset = c("all classified steps", "all classified", "all classified", "steady",
             "steady", "steady", "steady"),
  n = c(sum(stats::complete.cases(stepT[, PCA_FEATS_FIG])),
        sum(!is.na(stepT$gaitObjective)),
        sum(!is.na(stepT$gaitObjective)),
        sum(stepT$accClass == "steady", na.rm = TRUE),
        sum(stride$accClass == "steady" &
              as.character(stride$gaitObjective) %in% GAIT_LEVELS, na.rm = TRUE),
        sum(stepT$accClass == "steady", na.rm = TRUE),
        sum(stepT$accClass == "steady", na.rm = TRUE)))
readr::write_csv(fig_subsets, "output/figure_observation_subsets.csv")
cat("Observation subset behind each figure:\n"); print(fig_subsets)

readr::write_csv(main$tab, "output/GaitSummary_FourCell_Steady.csv")
readr::write_csv(main$n,   "output/GaitSummary_FourCell_Steady_cellN.csv")
readr::write_csv(sup,      "output/GaitSummary_FourCell_Unsteady.csv")
readr::write_csv(sup_n,    "output/GaitSummary_FourCell_Unsteady_cellN.csv")
readr::write_csv(comp,     "output/GaitSummary_FourCell_composition.csv")

# The aerial-pendular note quotes counts, so it is built from the composition table rather
# than written into the prose, where it would go stale the next time the pipeline runs.
PCA_FEATS_PROSE <- paste(unname(c(recovery = "pendular recovery", congruity = "KE-PE congruity",
                                  CoTmech = "mechanical cost of transport",
                                  dutyFactor = "duty factor")[PCA_FEATS]), collapse = ", ")
.aw  <- comp[comp$gait == unname(GAIT4_LABELS["pendularRun"]), ]
.oth <- comp[comp$gait != unname(GAIT4_LABELS["pendularRun"]), ]
AW_NOTE <- sprintf(paste0("The pendular run is the aerial + pendular cell, which the three-gait form ",
  "folds into aerial run; splitting it out makes aerial run %d rather than %d QC-passed steps. ",
  "Of its %d steps, %d%% are steady, against %d to %d%% steady in the three populated gaits, ",
  "so the birds do not hold it as a steady gait (see the composition sheet)."),
  sum(stepT$gait4 == "aerialRun", na.rm = TRUE),
  sum(stepT$hasFlight == 1, na.rm = TRUE),
  .aw$total, .aw$pct_steady, min(.oth$pct_steady), max(.oth$pct_steady))

NOTE_MAIN <- c(
  "Steady steps only: QC-passed steps whose parent stride is steady under the primary energy-grade criterion, |dE_CoM/(m g L)| <= 0.05. The first row gives the steady count and, after 'of', the count over all three steadiness classes, so the fraction of each cell that the summary rests on is explicit. The accelerating and decelerating remainder is summarized in the unsteady supplement table.",
  "Gait cells are the two qualitative classifier features crossed: whether an aerial phase occurs (the summed vertical force reaches zero) and whether CoM kinetic and gravitational potential energy fluctuate out of phase or in phase over the cycle (KE-PE congruity below or above 50% of the cycle). Each feature is a binary physical distinction split where the qualitative change occurs, so no threshold is tuned and the number of gaits follows from which cells are populated. The rotation-sense row is measured independently of these features and is reported as a result.",
  AW_NOTE,
  "Every cell is median (interquartile range, given as first quartile to third quartile). One format is used throughout because several of these measures are skewed or bounded (leg stiffness has a right tail from steps with near-zero estimated leg compression; duty factor is bounded near 0.5), so a median and IQR describe each cell without a mean being pulled by a tail.",
  "Decimal places are set from the data, not chosen per measure: each measure is rounded to the number of places needed to show the narrowest interquartile width among the four cells to two significant figures, taken over all steadiness classes so the steady and unsteady tables share one format. Reporting a median more finely than that would imply a precision the spread does not support.",
  "Collision angle is reported in DEGREES for readability. The theoretical equivalence between the force- and velocity-weighted collision angle and the dimensionless mechanical cost of transport (Lee et al. 2011, 2013) holds for the angle in RADIANS, which is the form the pipeline computes and the form used for the correlation quoted in the text; the degree values here are 180/pi times it.",
  "Touchdown mechanical energy is the total CoM mechanical energy at touchdown, dimensionless as E/(m g L0) = 0.5 u_TD^2 + h_TD/L0, with u_TD the sagittal CoM velocity magnitude at touchdown from the per-step cycle traces and h_TD the CoM height above the contacting foot from the virtual leg (legLen_TD sin alpha_TD). Steps whose reconstructed touchdown height falls outside 0.4 to 1.6 L0 are excluded, because foot-marker dropout and gap-fill contaminate legLen_TD in a small minority of steps; see output/touchdown_energy_check.csv and 03_step_descriptors.R.",
  sprintf("PC1 and PC2 are the CoM-dynamics principal components over %s, fitted across all classified steps with PC1 oriented positive at the walking end and PC2 positive for unsteady, identical to Table S1 and to the PC-space supplement figure.", PCA_FEATS_PROSE),
  "Leg stiffness is the dimensionless leg-spring stiffness k L0/(m g), equivalently body weights per leg length (Geyer, Seyfarth & Blickhan 2006), computed in MATLAB (GaitSelMulti_ExportTidyCSV/addLegStiffness) as k = Fmax/dL with the leg compression from the sinusoidal-GRF inversion of Blum, Lipfert & Seyfarth (2009) method C. Steps whose inverted compression is non-physical are excluded; see the cellN sheet for the n behind each cell.",
  "Leg stiffness is a model-bound descriptor: it is defined only relative to the spring-mass template it is derived from, and in the two grounded cells, where double support means both limbs share the net force, it is an effective two-limb stiffness rather than a single-limb property (McMahon 1985; Farley, Glasheen & McMahon 1993; Geyer, Seyfarth & Blickhan 2006). The walking column additionally rests on a small touchdown displacement and carries wider uncertainty than the aerial columns.",
  "For comparison, Birn-Jeffery et al. (2014) fitted a single constant k L0/(m g) = 15 for guinea fowl by work-optimal trajectory optimization on aerial running trials.")
NOTE_SUP <- c(
  "Same measures and conventions as the steady-step table, for the two unsteady classes of the same energy-grade criterion: accelerating (dE_CoM/(m g L) > 0.05) above, decelerating (< -0.05) below.",
  "Steps inherit the steadiness class of their parent stride, so an accelerating step is one taken within an accelerating stride.")

wb <- createWorkbook(); hdr <- createStyle(textDecoration = "bold", valign = "top")
addWorksheet(wb, "T_Steady");   writeData(wb, "T_Steady", main$tab, headerStyle = hdr)
writeData(wb, "T_Steady", data.frame(Notes = NOTE_MAIN), startRow = nrow(main$tab) + 3)
setColWidths(wb, "T_Steady", cols = 1:5, widths = c(46, 26, 26, 26, 26))
addWorksheet(wb, "TS_Unsteady"); writeData(wb, "TS_Unsteady", sup, headerStyle = hdr)
writeData(wb, "TS_Unsteady", data.frame(Notes = NOTE_SUP), startRow = nrow(sup) + 3)
setColWidths(wb, "TS_Unsteady", cols = 1:6, widths = c(20, 46, 26, 26, 26, 26))
addWorksheet(wb, "Composition");      writeData(wb, "Composition", comp, headerStyle = hdr)
addWorksheet(wb, "cellN_Steady");     writeData(wb, "cellN_Steady", main$n, headerStyle = hdr)
addWorksheet(wb, "cellN_Unsteady");   writeData(wb, "cellN_Unsteady", sup_n, headerStyle = hdr)
saveWorkbook(wb, "output/GaitSummary_FourCell.xlsx", overwrite = TRUE)

if (requireNamespace("knitr", quietly = TRUE)) {
  writeLines(c("# Four-cell CoM-dynamics summary tables", "",
               "Generated by _RAnalysis/R/20_table_gait4_summary.R.", "",
               "## Steady steps", "", knitr::kable(main$tab, format = "pipe"), "",
               "Notes:", "", paste0(seq_along(NOTE_MAIN), ". ", NOTE_MAIN), "",
               "## Unsteady steps: accelerating (upper), decelerating (lower)", "",
               knitr::kable(sup, format = "pipe"), "",
               "Notes:", "", paste0(seq_along(NOTE_SUP), ". ", NOTE_SUP), "",
               "## Composition of each cell by steadiness class", "",
               knitr::kable(comp, format = "pipe"), "",
               "## Sample size behind each cell, steady", "",
               knitr::kable(main$n, format = "pipe"), "",
               "## Sample size behind each cell, unsteady", "",
               knitr::kable(sup_n, format = "pipe"), ""),
             "output_internal/GaitSummary_FourCell.md")
}

# Leg stiffness against dimensionless speed, as INTERNAL evidence only.
#
# The manuscript reports leg stiffness descriptively, as median (IQR) per gait cell in Table 1, and
# makes no claim about how it changes with speed: it is a secondary measure, unrelated to the
# primary hypothesis, and the other descriptors carry no statistical test either. These statistics
# are kept because they are the evidence behind that decision, and because a linear R2 must not be
# read as flatness here (the relationship is not monotone, so a straight-line fit is near zero for
# a reason unrelated to constancy). They go to output_internal/, not output/.
r2 <- function(y, x) { ok <- is.finite(y) & is.finite(x)
  if (sum(ok) < 3) return(NA_real_); summary(lm(y[ok] ~ x[ok]))$r.squared }
kleg_speed <- bind_rows(
  data.frame(subset = "analysis sample", n = sum(is.finite(stepT$kLeg_n)),
             R2_kLeg_vs_u = round(r2(stepT$kLeg_n, stepT$meanSpeed_n), 3)),
  data.frame(subset = "steady steps", n = sum(is.finite(stepT$kLeg_n[stepT$accClass == "steady"])),
             R2_kLeg_vs_u = round(r2(stepT$kLeg_n[stepT$accClass == "steady"],
                                     stepT$meanSpeed_n[stepT$accClass == "steady"]), 3)),
  stepT %>% filter(!is.na(gait4)) %>% group_by(subset = as.character(gait4)) %>%
    summarise(n = sum(is.finite(kLeg_n)),
              R2_kLeg_vs_u = round(r2(kLeg_n, meanSpeed_n), 3), .groups = "drop"))
readr::write_csv(kleg_speed, "output_internal/kleg_vs_speed_r2.csv")

# The shape a linear R2 cannot see: rank correlation with speed overall, per gait and per
# individual, plus decile medians of stiffness and of its two components.
sp_rho <- function(y, x) { ok <- is.finite(y) & is.finite(x)
  if (sum(ok) < 10) return(c(rho = NA_real_, p = NA_real_, n = sum(ok)))
  ct <- suppressWarnings(cor.test(y[ok], x[ok], method = "spearman"))
  c(rho = unname(ct$estimate), p = ct$p.value, n = sum(ok)) }
kl <- stepT %>% filter(is.finite(kLeg_n), is.finite(meanSpeed_n))
rows <- list(data.frame(subset = "analysis sample", t(sp_rho(kl$kLeg_n, kl$meanSpeed_n))),
             data.frame(subset = "steady steps",
                        t(sp_rho(kl$kLeg_n[kl$accClass == "steady"],
                                 kl$meanSpeed_n[kl$accClass == "steady"]))))
for (g in c("walk", "groundedRun", "aerialRun")) {
  d <- kl %>% filter(as.character(gaitObjective) == g)
  rows[[length(rows) + 1]] <- data.frame(subset = g, t(sp_rho(d$kLeg_n, d$meanSpeed_n)))
}
kleg_rank <- bind_rows(rows) %>% mutate(rho = round(rho, 3), p = signif(p, 3))
readr::write_csv(kleg_rank, "output_internal/kleg_speed_rank.csv")
cat("Leg stiffness vs speed, Spearman rank correlation:\n"); print(as.data.frame(kleg_rank))

nb <- 10
kleg_bins <- kl %>%
  mutate(bin = dplyr::ntile(meanSpeed_n, nb)) %>%
  group_by(bin) %>%
  summarise(n = dplyr::n(), u_median = round(median(meanSpeed_n), 3),
            kLeg_median = round(median(kLeg_n), 2),
            kLeg_q1 = round(quantile(kLeg_n, .25), 2), kLeg_q3 = round(quantile(kLeg_n, .75), 2),
            peakF_BW_median = round(median(peakVertForce_BW, na.rm = TRUE), 3),
            legCompress_median = round(median(legCompress_n, na.rm = TRUE), 4), .groups = "drop")
readr::write_csv(kleg_bins, "output_internal/kleg_speed_bins.csv")

# Internal evidence, not a manuscript figure: leg stiffness against speed, with the decile median
# as the heavy line and the linear fit as the light line.
if (requireNamespace("ggplot2", quietly = TRUE)) {
  suppressPackageStartupMessages(library(ggplot2))
  klp <- kl %>% filter(as.character(gaitObjective) %in% GAIT_LEVELS) %>%
    mutate(gait = factor(as.character(gaitObjective), levels = GAIT_LEVELS))
  pk <- ggplot(klp, aes(meanSpeed_n, kLeg_n)) +
    geom_point(aes(colour = gait, shape = gait), alpha = 0.25, size = 0.9) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                colour = "grey55", linewidth = 0.6) +
    geom_line(data = kleg_bins, aes(u_median, kLeg_median),
              colour = "black", linewidth = 1.2, inherit.aes = FALSE) +
    geom_point(data = kleg_bins, aes(u_median, kLeg_median),
               colour = "black", size = 2, inherit.aes = FALSE) +
    scale_color_gait() + scale_shape_gait() +
    coord_cartesian(ylim = c(0, quantile(klp$kLeg_n, 0.99, na.rm = TRUE))) +
    labs(x = LAB_U, y = expression("leg stiffness  " * italic(k) * italic(L)[0]/(italic(m) * italic(g)))) +
    theme_daley()
  ggsave("output_internal/kleg_vs_speed.png", pk, width = 5.5, height = 3.6, dpi = 300, bg = "white")
  ggsave("output_internal/kleg_vs_speed.pdf", pk, width = 5.5, height = 3.6, device = grDevices::cairo_pdf)
  cat("wrote output_internal/kleg_vs_speed.{png,pdf}\n")
}
cat("Leg stiffness by speed decile:\n"); print(as.data.frame(kleg_bins))

# within-individual rank correlations, so the pooled value is not carried by between-bird spread
kleg_bird <- kl %>% group_by(subjectID) %>%
  filter(dplyr::n() >= 20) %>%
  summarise(n = dplyr::n(), rho = round(suppressWarnings(cor(kLeg_n, meanSpeed_n, method = "spearman")), 3),
            .groups = "drop")
readr::write_csv(kleg_bird, "output_internal/kleg_speed_rank_by_individual.csv")
cat(sprintf("within-individual rank correlation positive in %d of %d birds\n",
            sum(kleg_bird$rho > 0), nrow(kleg_bird)))

# gait differences, since "does not change substantially with gait" was also asserted
kw <- suppressWarnings(kruskal.test(kLeg_n ~ droplevels(factor(as.character(gaitObjective))),
                                    data = kl %>% filter(as.character(gaitObjective) %in%
                                      c("walk","groundedRun","aerialRun"))))
readr::write_csv(data.frame(test = "Kruskal-Wallis, kLeg_n by three-gait class",
                            statistic = round(unname(kw$statistic), 1), df = unname(kw$parameter),
                            p = signif(kw$p.value, 3)), "output_internal/kleg_gait_test.csv")
cat(sprintf("Kruskal-Wallis kLeg by gait: chi2 = %.1f, df = %d, p = %s\n",
            kw$statistic, kw$parameter, format.pval(kw$p.value)))
cat("Leg stiffness vs dimensionless speed (R2):\n"); print(as.data.frame(kleg_speed))

# What the leg-stiffness estimate is actually sensitive to. In the method-C inversion the
# normalised compression is (1 - sin(alpha_TD)) + C/L0, with C = (Fmax/m)(tc/pi)^2 - (g/8)tc^2
# independent of L0, so khat = (Fmax/BW) / that. The Discussion compares our values against two
# published estimates obtained by other conventions, and the elasticities below are what say how
# far such a comparison can be pushed. Emitted rather than asserted, and evaluated at each step
# then summarised, so the numbers describe this sample.
kleg_sens <- local({
  g <- 9.81
  d <- stepT %>% filter(is.finite(kLeg_n), is.finite(legAngle_TD), is.finite(L0_m),
                        is.finite(contactTime), is.finite(peakVertForce_BW))
  a  <- d$legAngle_TD * pi / 180
  Cc <- d$peakVertForce_BW * g * (d$contactTime / pi)^2 - (g / 8) * d$contactTime^2
  den <- (1 - sin(a)) + Cc / d$L0_m
  data.frame(
    quantity = c("touchdown leg angle alpha_TD (deg)",
                 "share of normalised compression from (1 - sin alpha_TD)",
                 "share of normalised compression from C/L0",
                 "d(ln kLeg)/d(ln L0)",
                 "d(ln kLeg)/d(alpha_TD), % per degree",
                 "steps used"),
    value = round(c(median(d$legAngle_TD),
                    median((1 - sin(a)) / den),
                    median((Cc / d$L0_m) / den),
                    median((Cc / d$L0_m) / den),
                    median(cos(a) / den) * pi / 180 * 100,
                    nrow(d)), 3))
})
readr::write_csv(kleg_sens, "output_internal/kleg_sensitivity.csv")

# Touchdown leg angle, over the ANALYSIS SAMPLE and with its spread, and split by cohort. Methods
# compares it against the 122.6 (5.4) deg Blum et al. 2014 report, and the offsets that set it are
# fitted per bird-session, so the cohort grain is the like-for-like one: the 2012 birds are the
# collections they measured. A pooled median that sits between two cohorts is not evidence of
# agreement with either.
alpha_tbl <- local({
  d <- stepT %>% filter(is.finite(legAngle_TD))
  one <- function(x, lab) data.frame(subset = lab, n = length(x),
                                     median_deg = round(median(x), 2), sd_deg = round(sd(x), 2),
                                     q1 = round(quantile(x, .25), 2), q3 = round(quantile(x, .75), 2))
  bind_rows(one(d$legAngle_TD, "analysis sample"),
            one(d$legAngle_TD[grepl("RVC", d$study)], "RVC 2008-2009 cohort"),
            one(d$legAngle_TD[!grepl("RVC", d$study)], "Blum 2012 cohorts"))
})
readr::write_csv(alpha_tbl, "output/touchdown_angle_by_cohort.csv")
cat("Touchdown leg angle (deg):\n"); print(as.data.frame(alpha_tbl))
cat("Leg-stiffness sensitivity:\n"); print(as.data.frame(kleg_sens))

# How often the spring-mass inversion returns a usable compression, by gait. It fails where the
# leg extends over stance, which is what vaulting does, so the loss is gait-dependent and the
# walking cell is summarised over a selected subset of its steps.
kleg_cov <- stepT %>% filter(!is.na(gaitObjective)) %>%
  group_by(gait = as.character(gaitObjective)) %>%
  summarise(steps = n(), with_kLeg = sum(is.finite(kLeg_n)),
            pct_with_kLeg = round(100 * with_kLeg / steps, 1), .groups = "drop")
readr::write_csv(kleg_cov, "output_internal/kleg_coverage_by_gait.csv")
cat("Leg-stiffness fit coverage by gait:\n"); print(as.data.frame(kleg_cov))

cat(sprintf("20_table_gait4_summary: four-cell tables (steady n = %d; accel n = %d; decel n = %d).\n",
            sum(stepT$accClass == "steady"), sum(stepT$accClass == "accelerating"),
            sum(stepT$accClass == "decelerating")))
cat("Leg stiffness by cell, steady steps (median (IQR), BW/L0):\n")
print(as.data.frame(stepT %>% filter(accClass == "steady") %>% group_by(gait = gait4) %>%
  summarise(n = sum(is.finite(kLeg_n)), kLeg = medIQR(kLeg_n, 1), .groups = "drop")))
cat("Cell composition by steadiness:\n"); print(as.data.frame(comp))
