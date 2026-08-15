# Guinea fowl gait selection: analysis code

Analysis code for a study of gait selection in guinea fowl (*Numida meleagris*) walking and
running overground on the level.

Center-of-mass motion is reconstructed from force-plate recordings by path-matched double
integration. Gait is classified from two binary mechanical features, whether an aerial phase
occurs and whether center-of-mass kinetic and potential energy fluctuate out of phase or in
phase. The rotation sense of the center-of-mass velocity loop, the hodograph, is then tested as a
criterion for the pendular-to-bouncing transition, against those features and against six
established gait descriptors.

**Paper:** Daley, M. A. and Birn-Jeffery, A. A velocity-loop view of avian gait: hodograph
rotation sense characterizes the walking-to-running continuum in guinea fowl. *Biology Open*
(submitted for review). The volume, article number and article DOI are added on acceptance.

**Data:** the recordings and the processed per-step and per-stride tables are archived on Dryad at
doi:10.5061/dryad.7m0cfxqc4. This repository is the code; the Dryad package is the data. Neither
contains the other, and the section below says how they meet.

**Version:** `v1.0.0`, cut from the working repository at commit
`6a8a77b72884ecf16528165500f15a25a9fabf74`. This repository keeps its revision history
from v1.0.0 forward: each release is tagged, and the commits between tags are the changes made
since the version the paper cites.

## Start here

`PIPELINE.md` is the order and output map: which script runs when, what each one writes, and
which element of the paper each output feeds. Read it before running anything.

Every number in the paper is emitted to a file by the pipeline. So a reader checking a value
looks for the output file that carries it, and `PIPELINE.md` says which script writes it.

## Layout

- `_MatlabProcessingCode/` the MATLAB phase: trial roster, import and filtering, whole-bout
  path-matched center-of-mass reconstruction, foot-contact gait-event detection, per-step and
  per-stride measures, the tidy-table export, the quality-control gate, and the exporters that
  build the Dryad package. `GaitSelMulti_RunFullWorkflow.m` is the master script and
  `helpers/` holds the biomechanics kernels, including `PathMatchedDoubleIntegration.m`.
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

The raw recordings are not in this repository. Download the Dryad package and place its
`per_trial_timeseries/` where `projectPaths.m` looks, or point that function at your copy: it
checks inside the project first and then one level above, so no path is hard-coded.

    matlab -batch "cd _MatlabProcessingCode; GaitSelMulti_RunFullWorkflow"

The master script runs the MATLAB phase and then hands off to R once. The quality-control gate is
computed in MATLAB so the Dryad export can embed it and the R phase reads rather than redefines
it, which is what keeps one definition of the analyzed set. Each `%%` block is independently
re-runnable, and an import cache means the raw recordings are read only when it is rebuilt
(`GaitSelMulti_BatchProcess(struct('useCache',false))`).

To run the R phase alone, from `_RAnalysis/`: `Rscript run_all.R`.

**One ordering note.** MATLAB STEP 5 embeds the gait and analysis-sample labels that R script
`04` writes, so on a first build run the sequence through and then re-run STEP 5, and then R
script `21` once more to restore the package's label columns. The full order is 1-5, 7, 5, 21.
`PIPELINE.md` gives the reason and the check that it has been done.

`_RAnalysis/R/test_kernels.R` holds known-answer tests for the computations a reader cannot check
by eye, including the shoelace signed area and its sign convention. It writes nothing and is not
part of `run_all.R`.

## The analysis sample

One gate, in `_RAnalysis/R/04_analysis_sample.R`: the MATLAB quality-control gate plus every
gait-space descriptor defined, giving **2588 steps and 945 strides** of the 3961 steps and 1731
strides detected, over 334 level trials from 30 recording sessions and 13 individuals. Nothing
downstream filters again, so any count reported as a fraction of all steps has 2588 as its
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
    version v1.0.0. https://github.com/NeuroMechLab/GuineaFowlGait-code

## Contact

Neuromechanics Lab (Daley), Department of Ecology and Evolutionary Biology, University of
California, Irvine. madaley@uci.edu
