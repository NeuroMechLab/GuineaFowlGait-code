# Pipeline order and output map

Every figure, table and number in the paper comes from the sequence below. Scripts are numbered
in the order they run, and the number is the order: MATLAB `STEP 1` to `STEP 7`, then R `01` to
`21`. Run the whole thing from
`_MatlabProcessingCode/GaitSelMulti_RunFullWorkflow.m`, which ends by calling the R phase.

## Requirements

MATLAB R2019b or later with the Statistics, Signal Processing and Curve Fitting toolboxes. R 4.6.0
with dplyr, ggplot2, patchwork, readr, mgcv, cluster, openxlsx and knitr. The raw force and marker
files in the `LEVEL_collation` archive are needed only to rebuild the import cache in STEP 2; see
Locating the raw archive below for where it is kept.

## Phase 0: the centre-of-mass offsets, in the shared core

This runs BEFORE the guinea fowl pipeline and lives outside this project, because the offsets are
fitted once for the whole archive and consumed by all three projects. It needs re-running only when
the roster, the mass records or the reconstruction change.

| Step | Script | Writes |
|:--|:---|:---|
| 0a | `AvianGaitCore/matlab/runFitComOffsets.m` | `LEVEL_collation/_metadata/CoM_offsets_fitted.csv`, `_PilotRun/comOffsetValidation/` |
| 0b | `AvianGaitCore/matlab/buildMorphologyOffsetTable.m` | `LEVEL_collation/_metadata/Morphology_offsets_lengths.csv` |

Step 0a fits the offset per bird-session from the net-zero pitch-impulse condition over whole strides, and
runs with the offsets disabled because it is fitting them. Method and rationale:
`_Documentation/PITCH_CORRECTION_METHOD.md`. Its validation figures need a person to look at them,
because the failure mode is anatomically implausible output rather than an error.

Add the core to the path with `coreOnPath(projectDir)`, never with a hand-ordered `addpath`. This
project vendors its own copy of five kernels, and whichever directory is added last wins.

## MATLAB phase: raw files to tidy tables

Run order is fixed by the dependency chain. The quality-control gate is computed here, once, so
that the data package can embed it and the R phase can read it rather than redefine it.

STEP 2 reads `CoM_offsets_fitted.csv` and applies the offset for that trial's bird-session where the
CoM proxy is built, so **the import cache must be rebuilt (`useCache=false`) after any change to the
offsets**. It also derives body mass per bird-session from the force record. STEP 3 derives L0 as
the median touchdown leg length over the bird-session's non-aerial steps, falling back to that
bird's sessions pooled and then to a within-sample allometry; `morphology.csv` records which branch
set each value.

Body mass, L0 and the CoM offset all key on one session key, `birdCode|YYYY-MM`. Month resolution
separates every guinea fowl session, since no bird was recorded on two dates within a month.

| Step | Script | Writes |
|:--|:---|:---|
| 1 | `GaitSelMulti_BuildRoster.m` | `_RAnalysis/data/trialRoster_multi.csv` (334 trials) |
| 2 | `GaitSelMulti_BatchProcess.m` | `SavedResults/perStepMeasures_multi.mat`, `perStep_raw_multi.csv`, `perStride_raw_multi.csv` |
| 3 | `GaitSelMulti_ExportTidyCSV.m` | `data/perStep_long.csv`, `perStride_long.csv`, `morphology.csv` |
| 4 | `GaitSelMulti_QCGate.m` | `data/step_qcpass.csv`, `stride_qcpass.csv`, `step_steadiness.csv` |
| 5 | `GaitSelMulti_ExportDryad.m` | Dryad `per_trial_timeseries/`, `trial_index.csv`, `steps_index.csv` |
| 5 | `GaitSelMulti_ExportCycleTraces.m` | `data/cycleTracesStep.csv` |
| 5 | `GaitSel_ExportDryadMeanCycles.m` | Dryad `mean_stepcycle_by_gait_speed/mean_cycles.{csv,mat}` |
| 6 | (QC diagnostic figures) | not part of this bundle; the workflow's `isfolder` guard skips the step |
| 7 | `_RAnalysis/run_all.R` | the R phase below |

Supporting kernels live in `_MatlabProcessingCode/helpers/`: `PathMatchedDoubleIntegration.m` (the
center-of-mass reconstruction), `scaleByDimension.m` (dynamic-similarity normalization),
`assembleForce.m` and `corrSafe.m` (fore-aft axis resolution), `filterForce.m`,
`detectStepsAccelMin.m`, `cleanFootTrack.m` (foot-marker spike removal),
`flagFreqLenOutliersSpeedLocal.m`, `flagSpeedDoubling.m` and `logMadFlag.m` (the outlier filters),
and `resolveRawPath.m` (raw-file location, below). `cleanFootSignal.m` is the single-channel form
of the same criterion and is not called by this pipeline.

Foot-marker tracks are spike-cleaned once, in `GaitSelMulti_ImportBout`, so the contact detector,
the virtual-leg geometry, the CoM-offset fit and the published per-trial series all measure the
same feet. Nothing downstream cleans them again.

### Locating the raw archive

`LEVEL_collation` is a shared multi-species archive, 3.6 GB over 9703 files, holding the raw
recordings for the guinea fowl collections and for the ostrich series the parallel projects use.
It is kept one level above this project so those projects read the same copy. `projectPaths`
searches for it under the project root and then one level above, accepting a candidate only if it
holds `OtherSpecies_byDate/Guinea fowl`, so either arrangement works. It reports the result as `P.collation`, empty when neither
confirms, since most of the pipeline runs from the tidy tables and never touches the raw files.
STEP 1 and the importer raise their own error naming both candidates.

The roster records every raw file by its path relative to `LEVEL_collation`, and `resolveRawPath`
rejoins it at read time to whichever location confirmed. Nothing stores an absolute path, so the
roster is portable and carries no home directory. `resolveRawPath` also accepts a project-relative
path and an absolute one, using the latter as given when it exists and otherwise re-resolving it
from its `LEVEL_collation` segment, which is what makes a roster built on another machine work
here.

STEP 5(a) reads `data/step_labels_analysis.csv`, which R script 04 writes, so `steps_index.csv`
carries the gait and analysis-sample columns only after STEP 5 has been re-run following the R
phase. On a first build, run the sequence through, then re-run STEP 5.

That second STEP 5 REVERTS two package files, so one more step closes the loop. STEP 5(a) copies
the tidy tables in as `perStep_long_multi.csv` and `perStride_long_multi.csv`, which drops the
label columns R script 21 adds, so script 21 runs again after it. The full order is 1-5, 7, 5, 21.
Script 21 is idempotent, building from `data/` rather than mutating the package copies, so the
repeat is safe. To check: `gait4` present in the package's `perStep_long_multi.csv` means the
relabel has run; absent means it has not. STEP 5(c) additionally reads
`output/speedbin_thresholds_R.csv` from the R phase and errors rather than guessing if it is
missing, so the package's mean cycles are binned on the same slow/fast threshold as the figures.

## R phase: tidy tables to figures, tables and statistics

`_RAnalysis/run_all.R` sources scripts 01 to 21 in this order. Script 22 is run separately, from
the project root, after the pipeline finishes. Files named `lib_*.R` hold shared definitions and
produce no output: `lib_theme.R` (house style and gait factor levels), `lib_format.R` (table
formatters), `lib_gam_crossing.R` (the single crossover-speed model), `lib_phase_register.R`
(Hilbert phase registration), `lib_hodograph.R` (closed-velocity-loop geometry: signed area, the
clockwise / counterclockwise area split, the retrograde turning fraction and the self-crossing
count, so that every loop descriptor reads the curve the classifier's signed area reads).

| # | Script | Manuscript element |
|:--|:---|:---|
| 01 | `01_load.R` | Methods 1 (individuals) |
| 02 | `02_clean.R` | Methods 3, 8 (quality control, reconstruction accuracy) |
| 03 | `03_step_descriptors.R` | Methods 7, Results 4 (collision angle, touchdown energy) |
| 04 | `04_analysis_sample.R` | Methods 2, Results 1, Table S2, Abstract (the analysis sample) |
| 05 | `05_mean_traces.R` | Results 7, Discussion 4 (vertical-force peak timing) |
| 06 | `06_fig_outlier_filter.R` | Fig. S1 |
| 07 | `07_steadiness.R` | Fig. S2, Methods 8, Results 1, 9 |
| 08 | `08_fig_hodograph_schematic.R` | Fig. 1 |
| 09 | `09_gait_continuity.R` | Tables S1, S3, Results 1, Discussion 2, Abstract |
| 10 | `10_fig_gaitspace_planes.R` | Fig. 2, Results 1 |
| 11 | `11_fig_energy_exchange.R` | Fig. 4, Methods 7, Results 4 |
| 12 | `12_hodograph_validation.R` | Results 3, Methods 11, Abstract (the three crossover speeds) |
| 13 | `13_fig_rotation_sense.R` | Fig. 3, Methods 6, Results 3 |
| 14 | `14_fig_hodographs.R` | Fig. 5, Results 5 |
| 15 | `15_sfig_trial_sequence.R` | Fig. S4, Results 5 |
| 16 | `16_fig_speed_relations.R` | Fig. 6, Results 6 |
| 17 | `17_fig_mean_traces.R` | Figs. 7, S5, S6 |
| 18 | `18_sfig_pc_steadiness.R` | Fig. S3, Table S1 |
| 19 | `19_table_gait_summary.R` | Results 2, 3 |
| 20 | `20_table_gait4_summary.R` | Tables 1, S4, Results 1 to 5, Methods 4, Discussion 5 |
| 21 | `21_dryad_labels.R` | the Dryad package's two tidy tables |

`run_all.R` lists each script's output files beside its `source()` call.

`R/test_kernels.R` is not sourced by `run_all.R` and writes nothing. It holds known-answer tests for
the computations a reader cannot check by eye: the shoelace signed area and its sign convention, the
method-C leg-compression inversion and its guards, the merged and split mis-cut flags, and the
crossover reader on monotone, non-monotone and exact-hit curves. Run it by hand
(`Rscript R/test_kernels.R`) after touching any of those; it exits non-zero on failure.

## Which figure comes from which script

| Figure | File | Script |
|:--|:---|:---|
| Fig. 1 | `Fig1_HodographSchematic.png` | 08 |
| Fig. 2 | `Fig2_GaitSpace.png` | 10 |
| Fig. 3 | `Fig3_RotationSense.png` | 13 |
| Fig. 4 | `Fig4_MetricMapping_Steady.png` | 11 |
| Fig. 5 | `Fig5_SteadyHodographs.png` | 14 |
| Fig. 6 | `Fig6_SpeedRelations_Steady.png` | 16 |
| Fig. 7 | `Fig7_MeanTraces_Steady.png` | 17 |
| Fig. S1 | `SFig_OutlierFilter.png` | 06 |
| Fig. S2 | `SFig_Steadiness.png` | 07 |
| Fig. S3 | `SFig_PCSteadiness.png` | 18 |
| Fig. S4 | `SFig_TrialSequence.png` | 15 |
| Fig. S5 | `SFig_MeanTraces_Accel.png` | 17 |
| Fig. S6 | `SFig_MeanTraces_Decel.png` | 17 |

The manuscript embeds the PNG; each script also writes a `cairo_pdf` companion of the same name
for submission. Figure numbers are not written into script names or into the figures themselves,
so a renumbering touches the captions and this table only.

| Table | Script | Artifact |
|:--|:---|:---|
| Table 1 | 20 | `GaitSummary_FourCell_Steady.csv`, `_Steady_cellN.csv` |
| Table S1 | 09 | `pca_summary.csv`, `pca_loadings_summary.csv` |
| Table S2 | 04 | `descriptor_correlations.csv` |
| Table S3 | 09 | `pam_clustering_summary.csv`, `pam_rotation_agreement.csv`, `_by_gait.csv` |
| Table S4 | 20 | `GaitSummary_FourCell_Unsteady.csv`, `_Unsteady_cellN.csv` |
| Table S5 | MATLAB STEP 3 | `data/morphology.csv` |

## Two output directories

`_RAnalysis/output/` holds only what the paper reports. Every file in it is either a figure the
manuscript embeds or an artifact behind a reported value. Seven
files reach the paper less directly, and are listed here so that nothing in this directory is
unaccounted for:

- `meanTraces_R.csv` and `meanTraces_annot_R.csv` are the plotted mean traces behind Figs. 7, S5
  and S6, and `speedbin_thresholds_R.csv` is the median-speed split those figures and Fig. 5
  share. They are figure source data rather than reported values.
- `GaitClassification_Tables.xlsx` and `GaitSummary_FourCell.xlsx` are cited in the numbers
  reference by sheet, as `xlsx T1` and `xlsx T_Steady`, not by filename.
- `GaitSummary_FourCell_Unsteady_cellN.csv` gives the per-measure denominator behind each cell of
  Table S4, which is what shows that only leg stiffness has an n below the cell count.
- `steadiness_by_gait.csv` is read by script 22 to build the per-gait steady counts it reports.

`_RAnalysis/output_internal/` holds products that are not reported. The leg-stiffness group is the
larger part of it: `kleg_speed_rank.csv`, `kleg_speed_bins.csv`,
`kleg_speed_rank_by_individual.csv`, `kleg_gait_test.csv`, `kleg_vs_speed_r2.csv` and
`kleg_vs_speed.{png,pdf}` are the evidence behind reporting leg stiffness descriptively and making
no claim about its speed dependence, since it is a secondary measure and the other descriptors
carry no statistical test either; `kleg_coverage_by_gait.csv` gives the per-gait share of steps
whose inverted compression is physical, and `kleg_sensitivity.csv` splits the normalized
compression between its touchdown-angle and leg-length terms, which is what shows the estimate is
driven by touchdown leg angle rather than by $L_0$. The rest are the within-bird correlation
spread, the steadiness-criterion comparison, the
clustering agreement against the classic criterion, the loop-area recomputation check, the
crossover-curve fits, the rotation-sense model note, `hodograph_crossover_check.png` (the loops
the self-crossing test flags and the loops it does not, drawn closed so the count can be checked
by eye), and a markdown mirror of the four-cell workbook. Each is kept as the evidence behind a
methods decision or a reviewer response. One
subdirectory sits alongside them, `dryad_readme_blocks/`, holding the three counted blocks that
the data package's README quotes, so that README can be rebuilt from the pipeline's own outputs
rather than hand-edited.

## Figure sizes

Every figure is exported at the width it is printed at, single column 88 mm or double column
183 mm, within a 210 mm height bound, as a vector PDF for the journal and a 300 dpi PNG for the
built manuscript. `bio_save()` in `lib_theme.R` writes both and refuses a canvas outside those
bounds; `BIO_W1`, `BIO_W2`, `BIO_HMAX`, `BIO_BASE` and `BIO_TAG` beside it are the specification,
and `theme_bio()` sets Arial at 8 pt with 12 pt bold capital panel tags.

Exporting at the printed width is what makes a nominal point size the printed point size. A canvas
wider than the placed width shrinks every label by the ratio between the two, so changing a
figure's width means changing the type size with it unless `bio_save()` is given the printed
width.

Two helpers carry the shared-axis layout. `bio_drop_x()` removes the tick labels, title and ticks
from a panel whose x axis is carried by the panel below it, which is what lets stacked rows sit
together; apply it only where the panels are drawn on the same x range, which for these figures is
set by a common `coord_cartesian(xlim = ...)`. `bio_inset_legend()` places a key inside the panel
in normalized coordinates.
