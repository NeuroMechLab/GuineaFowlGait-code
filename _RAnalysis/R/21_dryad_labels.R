# 21_dryad_labels.R — update the Dryad package's per-step and per-stride tables so they carry the
# gait labels the paper uses:
#   gait, gait4       the classification the paper uses (aerial phase x KE-PE congruity)
#   gaitHodo          the hodograph-based label, named for what it is, since the paper tests
#                     rotation sense against the classification
#   analysisSample    1 for the rows behind every number in the paper, 0 otherwise
#   exclusionReason   why a row is not analysed
# The MATLAB-written originals under _RAnalysis/data/ are untouched; only the package copies are
# rewritten, and only in their label columns.
suppressPackageStartupMessages({library(dplyr)})

DRY <- file.path("..", "GaitSel_DryadPackage_AllGF")
if (!dir.exists(DRY)) {
  cat("21_dryad_labels: no Dryad package directory; skipped.\n")
} else {
  lab <- readr::read_csv("data/step_labels_analysis.csv", show_col_types = FALSE)

  # This script REBUILDS the deposit from the MATLAB tidy tables, so it needs the raw ones, where
  # the hodograph label is still called gaitObjective. A reader reproducing the paper from the
  # published data has the DERIVED tables instead, in which that column is already gaitHodo and the
  # label columns are already present: dropping perStep_long_multi.csv in as data/perStep_long.csv
  # is the obvious thing to try, so this script detects that input rather than failing on the
  # missing gaitObjective column. Nothing in the paper depends on it, so it says so and skips.
  .raw <- names(readr::read_csv("data/perStep_long.csv", n_max = 0, show_col_types = FALSE))
  if (!"gaitObjective" %in% .raw) {
    cat(paste0("21_dryad_labels: data/perStep_long.csv is already the published table ",
               "(gaitObjective has been renamed gaitHodo), so the deposit needs no rebuild; ",
               "skipped. Every figure, table and statistic the paper reports comes from scripts ",
               "01 to 20, which have already run.\n"))
  } else {

  # Built FROM the MATLAB tidy tables in data/, not by mutating the package copies, so re-running
  # is idempotent: a second run reads the same source and produces the same output rather than
  # failing on already-renamed columns.
  ps <- readr::read_csv("data/perStep_long.csv", show_col_types = FALSE) %>%
    rename(gaitHodo = gaitObjective, mechHodo = mechClass) %>%
    left_join(lab, by = c("boutID", "stepIndex")) %>%
    mutate(analysisSample = ifelse(is.na(analysisSample), 0L, analysisSample),
           exclusionReason = ifelse(is.na(exclusionReason),
                                    "did not pass the QC gate", exclusionReason)) %>%
    relocate(gait, gait4, steadiness, analysisSample, exclusionReason, .after = stepIndex)
  readr::write_csv(ps, file.path(DRY, "perStep_long_multi.csv"))

  # The stride sample is the one 04_analysis_sample.R fixed: strides that pass the stride QC gate
  # AND whose two steps are both in the step sample. Deriving it from the step flags alone would
  # give a larger count, because that misses the stride-level gate, so the flag is taken from the
  # `stride` table the pipeline actually analysed.
  strideSamp <- stride %>% distinct(boutID, strideIndex) %>% mutate(inSamp = TRUE)
  st <- readr::read_csv("data/perStride_long.csv", show_col_types = FALSE) %>%
    rename(gaitHodo = gaitObjective, mechHodo = mechClass) %>%
    left_join(strideSamp, by = c("boutID", "strideIndex")) %>%
    mutate(analysisSample = as.integer(!is.na(inSamp))) %>% select(-inSamp) %>%
    relocate(analysisSample, .after = strideIndex)
  readr::write_csv(st, file.path(DRY, "perStride_long_multi.csv"))

  cat(sprintf(paste0("21_dryad_labels: per-step %d rows (%d analysed), per-stride %d rows ",
                     "(%d analysed); gait labels replaced with the classification the paper uses.\n"),
              nrow(ps), sum(ps$analysisSample), nrow(st), sum(st$analysisSample)))
  print(table(gait = ps$gait[ps$analysisSample == 1], useNA = "no"))
  }
}
