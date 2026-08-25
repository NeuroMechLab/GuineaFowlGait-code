# 01_load.R — read the tidy hand-off tables written by the MATLAB pipeline.
suppressPackageStartupMessages({library(readr); library(dplyr)})

DATA_DIR <- "data"
stride_raw <- readr::read_csv(file.path(DATA_DIR, "perStride_long.csv"), show_col_types = FALSE)
step_raw   <- readr::read_csv(file.path(DATA_DIR, "perStep_long.csv"),   show_col_types = FALSE)
morph      <- readr::read_csv(file.path(DATA_DIR, "morphology.csv"),     show_col_types = FALSE)

# Reconciled individual ID (random-effect grouping). The session-tagged `bird` is kept for the
# per-session physics (mass and L0 legitimately differ between sessions), but colour codes name
# the same individuals within a cohort, so the statistical random effect groups them: the RVC
# 2008-2009 collection is one cohort of birds, and the two 2012 Blum studies (Feb and June) are
# the same cohort. The bird the June 2012 records label `noc`, carrying no colour band, is the
# blue bird; the roster's individualCode carries the identification and the two
# records that establish it, and sets colourCode here.
add_subject <- function(df) {
  cohort <- ifelse(grepl("RVC", df$study), "rvc", "blum12")
  df$subjectID <- factor(paste0(cohort, "_", df$colourCode)); df
}
stride_raw <- add_subject(stride_raw); step_raw <- add_subject(step_raw)

cat(sprintf("Loaded %d strides, %d steps, %d bird-sessions, %d unique individuals.\n",
            nrow(stride_raw), nrow(step_raw), nrow(morph), nlevels(step_raw$subjectID)))
