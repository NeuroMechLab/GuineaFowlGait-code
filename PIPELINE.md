# Pipeline order and output map

Every figure, table and number in the paper comes from the sequence below, and the whole of it
runs from the Dryad data package. Scripts are numbered in the order they run: the MATLAB phase,
then R `01` to `21`.

## What you need

The data package, from Dryad at doi:10.5061/dryad.7m0cfxqc4. Unzip `per_trial_timeseries.zip` and
`mean_stepcycle_by_gait_speed.zip` inside it, then put the package at `GaitSel_DryadPackage_AllGF/`
beside `_MatlabProcessingCode/`. `projectPaths.m` locates it from its own file location, so no path
needs editing.

MATLAB R2019b or later with the Statistics, Signal Processing and Curve Fitting toolboxes. R 4.6.0
with dplyr, tidyr, ggplot2, patchwork, readr, mgcv, cluster, openxlsx and knitr; systemfonts is
used if present and skipped if not.

The raw force-plate and marker recordings are not part of this release and are not needed. The
deposited per-trial series carry the net ground reaction force, the time base, the kinematic
centre-of-mass proxy, both foot marker tracks, and body mass, leg length and g, which is
everything the measurement kernels read.

## MATLAB phase: the deposited trials to tidy tables

| Step | Script | Reads | Writes |
|:--|:---|:---|:---|
| 1 | `GaitSelMulti_BatchFromDryad.m` | the package's `per_trial_timeseries/` and `morphology_multi.csv` | `SavedResults/perStepMeasures_multi.mat`, `data/perStep_raw_multi.csv`, `data/perStride_raw_multi.csv`, `rerun/reconstruction_check.csv` |
| 2 | `GaitSelMulti_ExportTidyCSV.m` | `SavedResults/perStepMeasures_multi.mat` | `data/perStep_long.csv`, `perStride_long.csv`, `morphology.csv` (and the `_multi` names) |
| 3 | `GaitSelMulti_QCGate.m` | the tidy tables | `data/step_qcpass.csv`, `stride_qcpass.csv`, `step_steadiness.csv` |
| 4 | `GaitSelMulti_ExportCycleTraces.m` | the package's `per_trial_timeseries/` | `data/cycleTracesStep.csv` |

Step 1 rebuilds each trial into the bout struct the kernels expect and runs the same three of
them the paper used: `GaitSel_ReconstructCoM` for the path-matched double integration,
`GaitSel_DetectGaitEvents` with `GaitSel_FootContacts` to cut steps at foot-marker touchdowns, and
`GaitSel_PerStepStrideMeasures` for the per-step and per-stride reduction. `helpers/` holds the
kernels those three call, including `PathMatchedDoubleIntegration.m`.

Nothing writes into the data package. Step 1's outputs go to `rerun/` and `SavedResults/`, so a
second run starts from the archived state rather than from the first run's results.

`rerun/reconstruction_check.csv` compares the re-run against the values the package carries for
the same trial: the vertical centre-of-mass path, the drift, the work-energy identity ratio and the
number of steps cut. Over all 334 deposited trials it agrees to 0.000 mm on the CoM path and the
drift, to 0 on the work-energy ratio, and on the step count for every trial.

One column cannot be reproduced from the package. `dutyLimb`, the plate-based duty-factor
cross-check, is estimated from the per-plate vertical force, which is not deposited, so it comes
back empty. The duty factor the paper reports is the marker-based one and is unaffected.

The quality-control gate is computed in MATLAB, once, in step 3, so the R phase reads it rather
than redefining it. That is what keeps one definition of the analyzed set.

## Running it

    matlab -batch "addpath(genpath('_MatlabProcessingCode')); \
        GaitSelMulti_BatchFromDryad; GaitSelMulti_ExportTidyCSV; \
        GaitSelMulti_QCGate; GaitSelMulti_ExportCycleTraces"
    cd _RAnalysis && Rscript run_all.R

The R phase alone reproduces every figure, table and statistic once the four MATLAB steps have
written `_RAnalysis/data/`. Running the whole sequence from the package gives 3842 detected steps
and 1678 strides, 2620 steps and 1129 strides through the quality-control gate, and an analysis
sample of 2580 steps and 944 strides over 243 trials, which is what the paper reports.

`GaitSel_ExportDryadMeanCycles.m` rebuilds the package's own `mean_stepcycle_by_gait_speed/` from
the per-trial series. It reads `output/speedbin_thresholds_R.csv`, so it runs after the R phase,
and it is only needed to regenerate that part of the package.

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
