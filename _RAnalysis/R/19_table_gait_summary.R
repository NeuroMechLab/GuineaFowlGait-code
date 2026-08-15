# 19_table_gait_summary.R — gait-classification summary tables for the paper and for reviewers
# (GaitClassification_Tables.xlsx + CSVs): T1 per-gait descriptive statistics of the
# classifier features (aerial phase via duty factor; KE-PE congruity) and the traditional and
# proposed metrics.
#
# pct_CCW, the counterclockwise fraction of each gait, is a measured result here: gait is
# classified from the aerial phase and congruity (02_clean.R), so the rotation sense of the
# CoM velocity loop is independent of the label. The manuscript-facing version of the same
# breakdown, with the agreement statistics, is 12_hodograph_validation.R; this table repeats
# it for reviewers alongside the descriptor distributions.
suppressPackageStartupMessages({library(openxlsx); library(dplyr)})
if (!dir.exists("output")) dir.create("output")

fmt_msd <- function(x) sprintf("%.3f (%.3f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))

# T1: per-gait mean (s.d.) of the classifier and descriptor variables.
cand <- c(Froude = "Froude", dimensionless_speed_u = "meanSpeed_n", duty_factor = "dutyFactor",
          recovery_pct = "recovery", congruity_pct = "congruity", collision_angle = "collisionAngle",
          CoT_mech = "CoTmech", hodograph_area = "hodoArea_n")
cand <- cand[cand %in% names(step)]
gs <- step %>% filter(gaitObjective %in% GAIT_LEVELS) %>%
  mutate(gait = factor(gaitObjective, levels = GAIT_LEVELS)) %>% group_by(gait)
# pct_CCW is computed in its OWN summarise, so it reads the numeric hodoArea_n rather than the
# string column across(fmt_msd) creates under dplyr's sequential evaluation.
pctCCW <- gs %>% summarise(n_steps = n(),
                           pct_pendular_congruity = round(100 * mean(congruity < 50, na.rm = TRUE)),
                           pct_CCW = round(100 * mean(hodoArea_n > 0, na.rm = TRUE)), .groups = "drop")
msd <- gs %>% summarise(across(all_of(unname(cand)), fmt_msd), .groups = "drop")
gait_summary <- pctCCW %>% left_join(msd, by = "gait") %>%
  rename(!!!setNames(unname(cand), names(cand))) %>% arrange(gait)
readr::write_csv(gait_summary, "output/gait_classification_summary.csv")


wb <- createWorkbook()
addWorksheet(wb, "T1_GaitClassification"); writeData(wb, "T1_GaitClassification", gait_summary)
# S_PCA: principal-component loadings of the standardized CoM-dynamics features
# (written by 09_gait_continuity.R, which runs before this module); supports the continuum
# claim alongside the PAM silhouette profile in Table S3.
pca_path <- "output/pca_loadings_summary.csv"
if (file.exists(pca_path)) {
  pca_tbl <- readr::read_csv(pca_path, show_col_types = FALSE)
  addWorksheet(wb, "S_PCA"); writeData(wb, "S_PCA", pca_tbl)
}
# S_PAM: the silhouette profile (manuscript Table S3) and the agreement between the k = 2
# partition and hodograph rotation sense, both written by 09_gait_continuity.R.
for (s in list(c("S_PAM", "output/pam_clustering_summary.csv"),
               c("S_PAMagree", "output/pam_rotation_agreement.csv"))) {
  if (file.exists(s[2])) {
    addWorksheet(wb, s[1]); writeData(wb, s[1], readr::read_csv(s[2], show_col_types = FALSE))
  }
}
saveWorkbook(wb, "output/GaitClassification_Tables.xlsx", overwrite = TRUE)

cat("19_table_gait_summary: wrote GaitClassification_Tables.xlsx and CSVs.\n")
cat("Gait classification summary:\n"); print(as.data.frame(gait_summary))
