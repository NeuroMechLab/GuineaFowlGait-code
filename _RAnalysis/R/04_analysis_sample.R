# 04_analysis_sample.R — apply the final gate and fix THE analysis sample.
#
# One sample runs through every analysis in the project: the gait classification, the steadiness
# labels, the PCA and clustering, the transitions, and every table and figure. It is defined here,
# once, and nothing downstream filters on descriptor availability again.
#
# A step is in the sample when it passes the MATLAB QC gate (work-energy residual and
# duration-normalized drift, read in 02_clean.R) AND every continuous descriptor the paper reports
# as a gait-space axis is defined for it. The second condition matters because two of those
# descriptors depend on foot-marker quality that the force-based QC gate cannot see:
#   - the marker-based per-limb duty factor is undefined where foot contact could not be
#     determined from the markers;
#   - the touchdown energy is set aside where the reconstructed CoM height above the contacting
#     foot falls outside 0.4 to 1.6 L0, which flags foot-marker dropout or gap-fill;
#     one step fails both.
# Both are data-quality failures, so they belong in the gate. Gating here, rather than inside each
# script, is what gives every analysis one denominator.
#
# Leg stiffness is NOT in the gate. Its missingness is a model-fit failure, a non-physical inverted
# leg compression, rather than a defect in the measurement, and it is reported as a model-bound
# descriptor rather than as a gait-space axis. It keeps a per-measure count in the Table S1 notes.
#
# Strides. A stride is in the sample when it passes the stride QC gate and BOTH of its steps are in
# the step sample, so a stride is always composed of steps the paper analyses. Stride-grain
# analyses (the steady hodographs, the stride steadiness counts, the transitions) use that set.
#
# Output: output/analysis_sample.csv, the counts behind every denominator in the paper.
suppressPackageStartupMessages({library(dplyr)})
if (!dir.exists("output")) dir.create("output")

# the continuous descriptors that define the gait space (Table S3), all required
SAMPLE_DESCRIPTORS <- c("recovery", "congruity", "CoTmech", "collisionAngle", "E_TD_n", "dutyFactor")
# Reader-facing name for each descriptor, defined once here beside the set itself and used by
# 09_gait_continuity.R for the PCA table as well. The data-column names are for the code and the
# archived data; a table a reader reads gives the term, since these abbreviations appear nowhere
# else in the paper.
SAMPLE_DESCRIPTOR_LABELS <- c(recovery = "Pendular recovery", congruity = "KE-PE congruity",
                              CoTmech = "Mechanical cost of transport",
                              collisionAngle = "Collision angle", E_TD_n = "Touchdown energy",
                              dutyFactor = "Duty factor", hodoArea_signed = "Hodograph signed area")
stopifnot(all(SAMPLE_DESCRIPTORS %in% names(step)))

n_step_qc   <- nrow(step)
n_stride_qc <- nrow(stride)
keep_step <- stats::complete.cases(step[, SAMPLE_DESCRIPTORS])
lost_duty <- sum(!is.finite(step$dutyFactor) & is.finite(step$E_TD_n))
lost_etd  <- sum(is.finite(step$dutyFactor) & !is.finite(step$E_TD_n))
lost_both <- sum(!is.finite(step$dutyFactor) & !is.finite(step$E_TD_n))

# Authoritative per-step labels for the Dryad package, written BEFORE the sample is filtered so
# every detected step gets a row and the reason for any exclusion is explicit. The package's
# exporter reads this file for the gait label and the analysis flag, both of which live in R.
dryad_labels <- step %>%
  transmute(boutID, stepIndex,
            gait = as.character(gaitObjective),
            gait4 = as.character(gait4),
            steadiness = as.character(accClass),
            analysisSample = as.integer(keep_step),
            exclusionReason = dplyr::case_when(
              keep_step ~ "included in the analysis sample",
              !is.finite(dutyFactor) & !is.finite(E_TD_n) ~ "duty factor and touchdown height",
              !is.finite(dutyFactor) ~ "per-limb duty factor undetermined from foot markers",
              !is.finite(E_TD_n) ~ "touchdown height outside 0.4 to 1.6 L0",
              TRUE ~ "descriptor missing"))
readr::write_csv(dryad_labels, "data/step_labels_analysis.csv")
cat(sprintf("04_analysis_sample: wrote data/step_labels_analysis.csv (%d rows, %d in the sample).\n",
            nrow(dryad_labels), sum(dryad_labels$analysisSample)))

step <- step[keep_step, ]
step$s0 <- 2 * floor((step$stepIndex - 1) / 2) + 1

# strides whose two steps are both in the step sample
both_steps <- step %>% count(boutID, s0, name = "nsteps") %>% filter(nsteps == 2)
stride <- stride %>% semi_join(both_steps, by = c("boutID", "strideIndex" = "s0"))

# Both the disjoint counts and the two totals are emitted, because the Methods sentence needs the
# totals: a step failing both criteria is one exclusion, not two.
samp <- data.frame(
  quantity = c("steps passing the MATLAB QC gate",
               "steps dropped: marker-based duty factor undefined only",
               "steps dropped: touchdown height outside 0.4 to 1.6 L0 only",
               "steps dropped: both",
               "steps dropped: duty factor undefined, total",
               "steps dropped: touchdown height out of range, total",
               "steps dropped, total",
               "STEPS IN THE ANALYSIS SAMPLE",
               "trials contributing to the analysis sample",
               "strides passing the stride QC gate",
               "STRIDES IN THE ANALYSIS SAMPLE (both steps retained)",
               "steady steps", "accelerating steps", "decelerating steps",
               "steady strides"),
  value = c(n_step_qc, lost_duty, lost_etd, lost_both,
            lost_duty + lost_both, lost_etd + lost_both,
            lost_duty + lost_etd + lost_both, nrow(step),
            dplyr::n_distinct(step$boutID),
            n_stride_qc, nrow(stride),
            sum(step$accClass == "steady"), sum(step$accClass == "accelerating"),
            sum(step$accClass == "decelerating"),
            sum(stride$accClass == "steady")))
readr::write_csv(samp, "output/analysis_sample.csv")

cat(sprintf(paste0("04_analysis_sample: %d of %d QC-passed steps retained (%d lost: %d duty factor, ",
                   "%d touchdown height, %d both); %d of %d strides retained.\n"),
            nrow(step), n_step_qc, n_step_qc - nrow(step), lost_duty, lost_etd, lost_both,
            nrow(stride), n_stride_qc))
cat("Gait counts in the analysis sample (steps):\n"); print(table(step$gaitObjective))
print(table(step$gait4))
cat("Steadiness (steps):\n"); print(table(step$accClass))

# Correlation structure of the descriptor set. A dimensionality claim rests on how far these
# descriptors are already redundant, so the matrix is emitted rather than left implicit: the
# collision angle and the mechanical cost of transport are linked by an approximation (Lee et
# al.) and correlate accordingly, and recovery/congruity and touchdown-energy/duty-factor are
# each strongly related pairs. Written with fixed decimal places, so a column does not mix
# "1" and "-0.71" with "-0.125".
#
# The hodograph signed area is appended to the matrix but is NOT in SAMPLE_DESCRIPTORS, so it
# stays out of the PCA and the clustering: the loop is the measure tested against the structure
# those find. It is here so a reader can see how far the loop is already carried by the
# established descriptors.
step$hodoArea_signed <- step$hodoArea_n
CORR_VARS <- c(SAMPLE_DESCRIPTORS, "hodoArea_signed")
cm <- cor(step[, CORR_VARS])
# Emitted with the descriptors NAMED and the columns NUMBERED to the rows. A square matrix cannot
# carry seven full names across its header, and numbering the columns is the standard way round
# that: nothing is abbreviated, so the table needs no key. The emitted file takes the same form as
# the manuscript table, so the table is a transcription of it and can be checked against it.
cmlab <- paste(seq_along(CORR_VARS), unname(SAMPLE_DESCRIPTOR_LABELS[CORR_VARS]))
cmout <- data.frame(descriptor = cmlab, apply(cm, 2, function(z) sprintf("%.3f", z)),
                    check.names = FALSE, stringsAsFactors = FALSE)
names(cmout) <- c("Descriptor", as.character(seq_along(CORR_VARS)))
readr::write_csv(cmout, "output/descriptor_correlations.csv")
cat("Descriptor correlation matrix (output/descriptor_correlations.csv):\n"); print(round(cm, 3))

# Strongest correlation the signed area has with an established descriptor, positive or
# negative: the one value Results quotes from this matrix. The whole ranked row is emitted, so
# "strongest" can be checked rather than taken on trust.
r_hodo <- cm["hodoArea_signed", SAMPLE_DESCRIPTORS]
best <- data.frame(descriptor = names(r_hodo), pearson_r = round(unname(r_hodo), 3))
best <- best[order(-abs(best$pearson_r)), ]
best$rank_by_magnitude <- seq_len(nrow(best))
readr::write_csv(best, "output/hodoarea_descriptor_correlations.csv")
cat("Signed loop area against each descriptor, strongest first:\n"); print(best, row.names = FALSE)
