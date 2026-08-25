# run_all.R — GF gait-selection R pipeline orchestrator (single hand-off from MATLAB).
#
# Run from _RAnalysis/:  Rscript run_all.R
# Inputs, all written by the MATLAB phase:
#   data/perStep_long.csv, data/perStride_long.csv, data/morphology.csv   tidy measure tables
#   data/step_qcpass.csv, data/stride_qcpass.csv, data/step_steadiness.csv   the QC gate
#   data/cycleTracesStep.csv                                  per-step 100-point cycle traces
#
# The scripts are numbered in the order they run, and each one's outputs are listed beside it.
# output/ holds the figures, tables and statistics the paper reports. output_internal/ holds
# validation and methods-development products that the paper does not report; PIPELINE.md maps
# every output to the manuscript element that uses it.

setwd_if <- function() if (basename(getwd()) != "_RAnalysis" && dir.exists("_RAnalysis")) setwd("_RAnalysis")
setwd_if()
for (d in c("output", "output_internal")) if (!dir.exists(d)) dir.create(d)

# --- libraries: sourced for their definitions, no outputs -------------------------------------
source("R/lib_theme.R")            # theme_daley(), GAIT_LEVELS/LABELS/COLORS, axis-label constants
# lib_format.R, lib_gam_crossing.R, lib_phase_register.R and lib_hodograph.R are sourced by the
# scripts that need them (18, 11/12, 13 and 13/15 respectively), so each declares its own
# dependency. lib_hodograph.R holds the closed-velocity-loop geometry: signed area, the CW/CCW
# area split, the retrograde turning fraction and the self-crossing count.

# --- 01 to 05: build the analysis sample ------------------------------------------------------
source("R/01_load.R")              # stride_raw, step_raw, morph; the reconciled subjectID
source("R/02_clean.R")             # stride, step, with the MATLAB QC and steadiness labels
                                   #   output/reconstruction_drift_by_gait.csv
source("R/03_step_descriptors.R")  # collision angle in degrees, touchdown mechanical energy
                                   #   output/touchdown_energy_check.csv
source("R/04_analysis_sample.R")   # THE analysis sample, gated once; counts in analysis_sample.csv
                                   #   data/step_labels_analysis.csv  (read back by MATLAB STEP 5)
                                   #   output/analysis_sample.csv, descriptor_correlations.csv,
                                   #   output/hodoarea_descriptor_correlations.csv
source("R/05_mean_traces.R")       # group-mean step cycles from the cycle traces
                                   #   output/meanTraces_R.csv, meanTraces_annot_R.csv,
                                   #   output/speedbin_thresholds_R.csv, vertical_force_peak_timing.csv

# --- 06 to 07: sample-composition supplements -------------------------------------------------
source("R/06_fig_outlier_filter.R")# Fig. S1  output/SFig_OutlierFilter.{pdf,png}
source("R/07_steadiness.R")        # Fig. S2  output/SFig_Steadiness.{pdf,png}
                                   #   output/steadiness_{stride,step}_counts.csv, by_gait.csv,
                                   #   output/steadiness_speed_coverage.csv, _exclusion_by_speed.csv
                                   #   output_internal/steadiness_criteria_comparison.csv

# --- 08 to 18: the figure suite and the continuum analysis ------------------------------------
source("R/08_fig_hodograph_schematic.R") # Fig. 1  output/Fig1_HodographSchematic.{pdf,png}
source("R/09_gait_continuity.R")   # Tables S1 and S3: PCA and PAM clustering
                                   #   output/pca_summary.csv, pca_loadings_summary.csv,
                                   #   output/pam_clustering_summary.csv, pam_rotation_agreement.csv,
                                   #   output/pam_rotation_agreement_by_gait.csv, gaitspace_stats.csv
                                   #   output_internal/pam_classic_agreement.csv
source("R/10_fig_gaitspace_planes.R")# Fig. 2  output/Fig2_GaitSpace.{pdf,png}
                                   #   output/gaitspace_figure_stats.csv
source("R/11_fig_energy_exchange.R")# Fig. 4  output/Fig4_MetricMapping_Steady.{pdf,png}
                                   #   output/energy_exchange_stats.csv, collision_cot_relationship.csv,
                                   #   output/hodoarea_cot_relationship.csv, fig4_axis_clipping.csv
                                   #   output_internal/energy_exchange_within_bird.csv
source("R/12_hodograph_validation.R")# rotation sense tested against the classification
                                   #   output/hodo_validation_by_gait.csv, hodo_validation_boundary.csv
                                   #   output_internal/hodo_validation_boundary_fits.csv
source("R/13_fig_rotation_sense.R")# Fig. 3  output/Fig3_RotationSense.{pdf,png}
                                   #   output/hodograph_reversal_diagnostics.csv, _overall.csv,
                                   #   output/hodograph_rotation_by_gait.csv, rotation_reversal.csv
                                   #   output_internal/hodoArea_recomputation_check.csv,
                                   #   output_internal/hodograph_rotation_model_note.txt
source("R/14_fig_hodographs.R")    # Fig. 5  output/Fig5_SteadyHodographs.{pdf,png}
                                   #   output/fig4_rotation_sense_observed.csv, fig5_loop_closure.csv
source("R/15_sfig_trial_sequence.R")# Fig. S4  output/SFig_TrialSequence.{pdf,png}
                                   #   output/trial_sequence_candidates.csv,
                                   #   output/trial_sequence_strides.csv
source("R/16_fig_speed_relations.R")# Fig. 6  output/Fig6_SpeedRelations_Steady.{pdf,png}
                                   #   output/gait_speed_overlap_steady.csv, _histogram_steady.csv
source("R/17_fig_mean_traces.R")   # Fig. 7 and Figs. S5, S6
                                   #   output/Fig7_MeanTraces_Steady.{pdf,png}
                                   #   output/SFig_MeanTraces_{Accel,Decel}.{pdf,png}
source("R/18_sfig_pc_steadiness.R")# Fig. S3  output/SFig_PCSteadiness.{pdf,png}

# --- 19 to 21: tables and the data-package labels ---------------------------------------------
source("R/19_table_gait_summary.R")# output/gait_classification_summary.csv,
                                   #   output/GaitClassification_Tables.xlsx
source("R/20_table_gait4_summary.R")# Tables 1 and S4, and the leg-stiffness rows
                                   #   output/GaitSummary_FourCell_{Steady,Unsteady}.csv and _cellN,
                                   #   output/GaitSummary_FourCell_composition.csv, .xlsx,
                                   #   output/figure_observation_subsets.csv,
                                   #   output/touchdown_angle_by_cohort.csv
                                   #   output_internal/GaitSummary_FourCell.md, kleg_*.csv,
                                   #   output_internal/kleg_vs_speed.{png,pdf}
source("R/21_dryad_labels.R")      # writes the Dryad package's two tidy tables
                                   #   ../GaitSel_DryadPackage_AllGF/perStep_long_multi.csv,
                                   #   ../GaitSel_DryadPackage_AllGF/perStride_long_multi.csv

cat("\nR pipeline complete. Paper outputs in _RAnalysis/output/,",
    "validation products in _RAnalysis/output_internal/.\n")
