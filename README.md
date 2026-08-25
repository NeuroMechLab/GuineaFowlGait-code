# Guinea fowl gait selection: analysis code

Analysis code for a study of gait selection in guinea fowl (*Numida meleagris*) walking and
running overground on the level.

Center-of-mass motion is reconstructed from force-plate recordings by path-matched double
integration. Gait is classified from two binary mechanical features, whether an aerial phase
occurs and whether center-of-mass kinetic and potential energy fluctuate out of phase or in
phase. The rotation sense of the center-of-mass velocity loop, the hodograph, is then tested as a
criterion for the pendular-to-bouncing transition, against those features and against six
established gait descriptors.

**Paper:** Daley, M. A. and Birn-Jeffery, A. (2026). A velocity-loop view of avian gait:
hodograph rotation sense characterizes the walking-to-running continuum in guinea fowl.
*Biology Open*, doi:10.1242/bio.062880

**Data:** the recordings and the processed per-step and per-stride tables are archived on Dryad at
doi:10.5061/dryad.7m0cfxqc4. This repository is the code; the Dryad package is the data. Neither
contains the other, and the section below says how they meet.

**Version:** `v1.1.0`, cut from the working repository at commit
`65492fb228615a64f19fd1aa17a5bc011aceedc8`. This repository keeps its revision history
from v1.0.0 forward: each release is tagged, and the commits between tags are the changes made
since the version the paper cites.

## Start here

`PIPELINE.md` is the order and output map: which script runs when, what each one writes, and
which element of the paper each output feeds. Read it before running anything.

Every number in the paper is emitted to a file by the pipeline. So a reader checking a value
looks for the output file that carries it, and `PIPELINE.md` says which script writes it.

## Layout

- `_MatlabProcessingCode/` the MATLAB phase, which runs from the Dryad package: whole-bout
  path-matched center-of-mass reconstruction, foot-contact gait-event detection, per-step and
  per-stride measures, the tidy-table export and the quality-control gate.
  `GaitSelMulti_BatchFromDryad.m` is the entry point and `helpers/` holds the biomechanics
  kernels, including `PathMatchedDoubleIntegration.m`.
- `_RAnalysis/` the R phase, orchestrated by `run_all.R`: scripts `01` to `21`, numbered in the
  order they run. `lib_*.R` files hold shared definitions. Inputs in `data/`, the figures, tables and statistics the paper reports in
  `output/`, and validation products it does not report in `output_internal/`.
- `GaitSel_DryadPackage_AllGF/README.md` the data package's own documentation: every file, every
  column, the units and the selection criteria. The data itself is on Dryad, not here.

This repository is the analysis code. The manuscript's own scaffolding lives with the manuscript
and is not here: the source, the bibliography and citation style, the build and prose-checking
scripts, the journal formatting filters, and the per-number provenance table. None of it bears on
how the science was computed, which is what `PIPELINE.md` and the two script directories set out.

## Requirements

MATLAB R2019b or later, with the Statistics, Signal Processing and Curve Fitting toolboxes. R
4.6.0, with `dplyr`, `tidyr`, `ggplot2`, `patchwork`, `readr`, `mgcv`, `cluster`, `openxlsx` and
`knitr`; `systemfonts` is used if present and skipped if not, to pick up the figures' typeface.
Building the manuscript additionally needs pandoc 3 or later, and a TeX installation with
`xelatex` for the PDF. The pipeline as reported was run on macOS with MATLAB R2025b and R 4.6.0.

`PIPELINE.md` is where this list is maintained; if the two ever disagree, that file is the one
kept current.

## Running it

Download the data package from Dryad, doi:10.5061/dryad.7m0cfxqc4. Unzip
`per_trial_timeseries.zip` and `mean_stepcycle_by_gait_speed.zip` inside it, then put the package
at `GaitSel_DryadPackage_AllGF/` beside `_MatlabProcessingCode/`. `projectPaths.m` finds it from
its own location, so nothing needs editing.

    matlab -batch "addpath(genpath('_MatlabProcessingCode')); \
        GaitSelMulti_BatchFromDryad; GaitSelMulti_ExportTidyCSV; \
        GaitSelMulti_QCGate; GaitSelMulti_ExportCycleTraces"
    cd _RAnalysis && Rscript run_all.R

The MATLAB phase re-runs the reconstruction, the gait-event detection and the per-step reduction
from the deposited per-trial series, then writes the tidy tables and the quality-control gate the
R phase reads. The gate is computed once, in MATLAB, so R reads rather than redefines it, which is
what keeps one definition of the analyzed set.

Nothing writes into the data package. The re-run's outputs go to `rerun/` and
`_MatlabProcessingCode/SavedResults/`, so a second run starts from the archived state.

`rerun/reconstruction_check.csv` compares the re-run against the values the package carries for
each trial. Over all 334 deposited trials it agrees to 0.000 mm on the vertical center-of-mass
path and on the drift, to 0 on the work-energy identity ratio, and on the number of steps cut for
every trial.

The raw force-plate and marker recordings are not part of this release and are not needed. One
column depends on them: `dutyLimb`, the plate-based duty-factor cross-check, comes back empty. The
duty factor the paper reports is the marker-based one and is unaffected.

To run the R phase alone, once `_RAnalysis/data/` has been written: `cd _RAnalysis && Rscript run_all.R`.

`_RAnalysis/R/test_kernels.R` holds known-answer tests for the computations a reader cannot check
by eye, including the shoelace signed area and its sign convention. It writes nothing and is not
part of `run_all.R`.

## The analysis sample

One gate, in `_RAnalysis/R/04_analysis_sample.R`: the MATLAB quality-control gate plus every
gait-space descriptor defined, giving **2580 steps and 944 strides** of the 3842 steps and 1678
strides detected, over 323 level trials from 30 recording sessions and 13 individuals. Nothing
downstream filters again, so any count reported as a fraction of all steps has 2580 as its
denominator. The Dryad package carries an `analysisSample` flag that reproduces exactly that set,
which is how a reader reproduces the paper's sample without rerunning the gate.

## Conventions

Any three-vector is ordered [medio-lateral, fore-aft (direction of travel), vertical], and units
are SI. Speed is reported dimensionless as *u* = *v*/sqrt(*gL*<sub>0</sub>), with Froude number
*u*<sup>2</sup>, where *L*<sub>0</sub> is the bird's leg length. Body mass and *L*<sub>0</sub> are
kept per recording session, because a bird's weight varied between sessions.

## License

Code is released under the MIT license, in `LICENSE`.
The Dryad package carries its own license terms for the data.

## Citing this

Cite the paper. If you need to cite the code specifically, cite the archived release:

    Daley, M. A. and Birn-Jeffery, A. (2026). Guinea fowl gait selection: analysis code,
    version v1.1.0. https://github.com/NeuroMechLab/GuineaFowlGait-code

## Contact

Neuromechanics Lab (Daley), Department of Ecology and Evolutionary Biology, University of
California, Irvine. madaley@uci.edu
