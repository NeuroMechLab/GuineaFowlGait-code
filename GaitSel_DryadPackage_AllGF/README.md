## Guinea fowl center-of-mass dynamics across the walking-to-running speed range

Data and code associated with the paper: **"A velocity-loop view of avian gait: hodograph rotation sense characterizes the walking-to-running continuum in guinea fowl"** Authors: Monica A. Daley and Aleksandra Birn-Jeffery\
*Biology Open* (2026). doi:10.1242/bio.062880

\***Author for correspondence:** [madaley@uci.edu](mailto:madaley@uci.edu)

## Version

This version exchanges the right and left foot-marker tracks in 41 trials whose marker labels were swapped in the original captures, and corrects this README's description of the coordinate frame and of several columns; no step, stride, measure or reported value changes. Version 2 corrected the sign of the horizontal ground reaction force in 77 trials and every quantity computed from it. Version History, below, describes both, and the effect of version 2 on the results the article reports, which were computed from the first version.

## Overview

This package contains ground reaction force recordings and reconstructed center-of-mass (CoM) trajectories for guinea fowl (*Numida meleagris*) running over level ground, compiled from three Royal Veterinary College data collections and reprocessed through one common pipeline. The analysis asks how the rotation sense of the CoM velocity loop, the hodograph, relates to the transition from pendular (walking) to bouncing (running) dynamics, and compares it against established center-of-mass criteria across a speed range that includes extensive grounded running.

Each of 334 trials is provided as the time series of its on-plate bout, the span in which the bird is both tracked and on the plates, in a standardized, bird-centered coordinate frame, together with tidy per-step and per-stride measure tables, per-step gait classification, and ensemble mean step cycles by gait and speed.

## Directory Structure

```
GaitSel_DryadPackage_AllGF/
├── README.md                                      # This file
├── trial_index.csv                                # One row per trial (334)
├── steps_index.csv                                # One row per step (3842)
├── perStep_long_multi.csv                         # Tidy per-step measures (3842 x 90)
├── perStride_long_multi.csv                       # Tidy per-stride measures (1678 x 78)
├── morphology_multi.csv                           # Per-session body mass and leg length (30)
├── trialRoster_multi.csv                          # Source file provenance per trial (334)
├── force_sign_corrections.csv                     # Force columns whose sign version 2 corrected (78)
├── foot_side_corrections.csv                      # Trials flagged by the foot-marker side check (69)
├── per_trial_timeseries.zip                       # Whole-bout time series, one pair per trial
│   └── per_trial_timeseries/                      #   unzips to this folder
│       ├── [bird]_[YYYYMMDD]_[trialStem].csv      # 20 channels at 500 Hz
│       └── [bird]_[YYYYMMDD]_[trialStem].mat      # Same, plus metadata and a step table
└── mean_stepcycle_by_gait_speed.zip               # Ensemble mean step cycles
    └── mean_stepcycle_by_gait_speed/              #   unzips to this folder
        ├── mean_cycles.csv                        # Long format, 10800 rows
        └── mean_cycles.mat                        # Same content as a MATLAB struct
```

The two folders are deposited as `.zip` archives because Dryad stores a flat file list and would
otherwise drop the folder they belong to. Each archive unzips to the folder named above, so the
layout after unzipping is the one the pipeline and the paths below expect. Everything else is a
flat file at the top level.

## Data Files

### File Naming Convention

Per-trial files follow the naming pattern:

```
[bird]_[Date]_[trialStem].csv
```

**Example:** `rvc_bly_20090505_L00cm_0009.csv`

| Component | Format         | Description                                                        |
| :-------- | :------------- | :----------------------------------------------------------------- |
| bird      | study_colour   | Study-session bird code (e.g. `rvc_bly`, `pot_red`, `surf_noc`)    |
| Date      | YYYYMMDD       | Recording date (e.g. 20090505 = 5 May 2009)                        |
| trialStem | L00cm_#### or LE0cm_#### | Terrain code and sequential trial number within that session |

The same string is the `boutID`, the join key across the per-trial, trial, step and stride files; `morphology_multi.csv` joins on `session`.

### Studies Included

All level (unperturbed) running trials with raw force and synchronized marker capture.

| Study tag                    | Trials | Steps | Strides | Dates                    | Plates | Description                                                     |
| :--------------------------- | -----: | ----: | ------: | :----------------------- | :----- | :-------------------------------------------------------------- |
| `RVC_Daley_2008-2011` | 189 | 1969 | 850 | 2008-11-20 to 2009-11-19 | 5 or 6 | Gait-selection collection plus the level baselines of an obstacle and a drop-perturbation collection |
| `Blum_DropVsPothole_2012-02` | 77 | 1049 | 466 | 2012-02-14 to 2012-02-17 | 6 | Level running, the control condition of a drop-versus-pothole study |
| `Blum_Surface_2012-06` | 68 | 824 | 362 | 2012-06-25 to 2012-06-28 | 6 | Level running on bare force plates, the flat rigid reference of a surface-compliance study |
| **Total** | 334 | 3842 | 1678 | | | |

Force was sampled at 500 Hz.

### Subject Information

Thirty bird-sessions over thirteen individuals. Colour codes name the same individual within a collection cohort but were reused across cohorts, so `subjectID` is the reconciled identity and is the correct grouping factor for any analysis pooling across studies. The two 2012 Blum studies, four months apart, are one cohort. The June bird recorded as `noc`, carrying no colour band, is the blue bird: the study's weight record lists blue among the seven birds run and the archive holds no `blu` for that study, it is the only unmatched bird, and the original analysis's published data uses blue's weighed mass for 'noc'. Body mass and leg length are set per session. Body mass is estimated from force: the median, over the session's trials, of the mean summed vertical ground reaction force in the analysis window, divided by g. The masses in the collection metadata, one value per bird in each study, are `bodyMass_kg` in `trialRoster_multi.csv`.

| Session | subjectID | Study | Session date | Mass (kg) | L0 (m) | L0 source | Trials | Steps | Strides |
| :------ | :-------- | :---- | :----------- | --------: | -----: | :-------- | -----: | ----: | ------: |
| `pot_bla` | `blum12_bla` | Blum_DropVsPothole_2012-02 | 2012-02 | 1.158 | 0.208 | session | 10 | 169 | 77 |
| `pot_blu` | `blum12_blu` | Blum_DropVsPothole_2012-02 | 2012-02 | 1.087 | 0.205 | session | 10 | 133 | 60 |
| `pot_g00` | `blum12_g00` | Blum_DropVsPothole_2012-02 | 2012-02 | 1.313 | 0.207 | session | 10 | 119 | 52 |
| `pot_gry` | `blum12_gry` | Blum_DropVsPothole_2012-02 | 2012-02 | 1.377 | 0.222 | session | 10 | 140 | 61 |
| `pot_reb` | `blum12_reb` | Blum_DropVsPothole_2012-02 | 2012-02 | 1.177 | 0.219 | session | 10 | 134 | 60 |
| `pot_red` | `blum12_red` | Blum_DropVsPothole_2012-02 | 2012-02 | 1.243 | 0.204 | session | 10 | 181 | 83 |
| `pot_rre` | `blum12_rre` | Blum_DropVsPothole_2012-02 | 2012-02 | 1.253 | 0.204 | allometry | 7 | 52 | 20 |
| `pot_yel` | `blum12_yel` | Blum_DropVsPothole_2012-02 | 2012-02 | 1.172 | 0.198 | session | 10 | 121 | 53 |
| `surf_bla` | `blum12_bla` | Blum_Surface_2012-06 | 2012-06 | 1.453 | 0.222 | session | 10 | 190 | 87 |
| `surf_g00` | `blum12_g00` | Blum_Surface_2012-06 | 2012-06 | 1.740 | 0.223 | session | 9 | 116 | 52 |
| `surf_gry` | `blum12_gry` | Blum_Surface_2012-06 | 2012-06 | 1.492 | 0.232 | session | 10 | 111 | 48 |
| `surf_noc` | `blum12_blu` | Blum_Surface_2012-06 | 2012-06 | 1.361 | 0.210 | session | 10 | 117 | 52 |
| `surf_red` | `blum12_red` | Blum_Surface_2012-06 | 2012-06 | 1.332 | 0.214 | session | 10 | 129 | 58 |
| `surf_rre` | `blum12_rre` | Blum_Surface_2012-06 | 2012-06 | 1.424 | 0.206 | allometry | 9 | 68 | 26 |
| `surf_yel` | `blum12_yel` | Blum_Surface_2012-06 | 2012-06 | 1.302 | 0.227 | session | 10 | 93 | 39 |
| `rvc_blu` | `rvc_blu` | RVC_Daley_2008-2011 | 2008-11 | 1.161 | 0.187 | session | 2 | 18 | 8 |
| `rvc_blu` | `rvc_blu` | RVC_Daley_2008-2011 | 2009-05 | 1.378 | 0.200 | session | 13 | 113 | 48 |
| `rvc_blu` | `rvc_blu` | RVC_Daley_2008-2011 | 2009-11 | 1.286 | 0.199 | bird | 10 | 42 | 14 |
| `rvc_bly` | `rvc_bly` | RVC_Daley_2008-2011 | 2008-11 | 1.277 | 0.197 | bird | 4 | 14 | 4 |
| `rvc_bly` | `rvc_bly` | RVC_Daley_2008-2011 | 2009-05 | 1.653 | 0.197 | session | 29 | 471 | 214 |
| `rvc_bly` | `rvc_bly` | RVC_Daley_2008-2011 | 2009-11 | 1.273 | 0.197 | bird | 11 | 80 | 32 |
| `rvc_red` | `rvc_red` | RVC_Daley_2008-2011 | 2008-11 | 1.380 | 0.211 | session | 4 | 29 | 12 |
| `rvc_red` | `rvc_red` | RVC_Daley_2008-2011 | 2009-05 | 1.934 | 0.197 | session | 20 | 286 | 128 |
| `rvc_red` | `rvc_red` | RVC_Daley_2008-2011 | 2009-11 | 1.338 | 0.199 | bird | 12 | 77 | 30 |
| `rvc_rey` | `rvc_rey` | RVC_Daley_2008-2011 | 2008-11 | 1.317 | 0.192 | session | 1 | 9 | 4 |
| `rvc_rey` | `rvc_rey` | RVC_Daley_2008-2011 | 2009-05 | 1.771 | 0.192 | session | 24 | 363 | 161 |
| `rvc_rey` | `rvc_rey` | RVC_Daley_2008-2011 | 2009-11 | 1.276 | 0.208 | session | 15 | 107 | 44 |
| `rvc_yel` | `rvc_yel` | RVC_Daley_2008-2011 | 2008-11 | 1.203 | 0.163 | session | 4 | 31 | 13 |
| `rvc_yel` | `rvc_yel` | RVC_Daley_2008-2011 | 2009-04 | 1.504 | 0.197 | session | 29 | 284 | 122 |
| `rvc_yel` | `rvc_yel` | RVC_Daley_2008-2011 | 2009-11 | 1.190 | 0.196 | bird | 11 | 45 | 16 |

Leg length `L0` is the median CoM-to-foot distance at touchdown over the session's steps with no aerial phase and a touchdown leg length between 0.10 and 0.40 m. Where a session yielded fewer than five such steps it falls back to that bird's sessions pooled, and for the one bird that never yielded five, to a within-sample allometry on body mass; the `L0 source` column records which applies.

### Data Selection

Every detected step, and every stride the pipeline forms, is provided, with flags rather than deletion, so a reader can reproduce the paper's sample or define a different one.

```
3842 steps detected
2618 pass the quality-control gate           (qcPass = 1)
   0 duty factor undetermined
  40 touchdown height out of range
   0 both
2578 in the analysis sample                  (analysisSample = 1)
```

| Flag             | Meaning                                                                                                                                                                                                                                     |
| :--------------- | :------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `outlier`        | 1 marks a mis-cut cycle: a log-MAD filter on step frequency and length judged within dimensionless-speed neighborhoods, combined with a speed-conditional filter catching merged (half-frequency) and split (double-frequency) cycles       |
| `qcPass`         | 1 marks steps passing the retention gate: `outlier = 0`, finite positive speed, work-energy identity residual at or below 0.20, and vertical reconstruction drift within a duration-normalized bound (the larger of 45 mm and 19.5 mm/s times the summed duration of the trial's steps, which is 45 mm at the median trial) |
| `analysisSample` | 1 marks the rows behind every number in the paper. Adds to the retention gate the requirement that every continuous gait-space descriptor be defined                                                                                        |
| `exclusionReason`| Names the reason for every row outside the analysis sample                                                                                                                                                                                  |

Two gait-space descriptors depend on foot-marker quality that cannot be detected based on a force-based gate: the per-limb duty factor is undefined where foot contact could not be determined from the markers, and the touchdown energy is considered infeasible where the CoM height above the contacting foot at touchdown, h = `legLen_TD_n` × sin(`legAngle_TD`) in leg lengths, is not strictly between 0.4 and 1.6.

**`analysisSample = 1` is the column to filter on to reproduce the paper.** It selects 2578 steps and 942 strides. Within the step sample the steadiness split is 929 steady, 942 accelerating and 707 decelerating. Of 1678 strides detected, 1126 pass the stride gate and 942 are in the analysis sample, those being the strides that pass the gate and whose two steps are both in the step sample.

Strides tile each trial's steps in pairs from the first step: the stride with `strideIndex` k spans steps k and k + 1, for odd k. A trial with an odd number of steps leaves its last step unpaired (160 trials), and the pipeline forms no stride from the final pair of steps in a trial with an even number of steps (163 trials, 326 steps), so 486 steps belong to no stride; their steadiness is graded on the step itself (see `steadiness`).

---

## Missing Values

Every CSV in this package writes a missing value as `NA`, and no cell is left blank. A value is
`NA` where the quantity is not defined for that row, and the flag or label that says why it is
undefined sits in the same row. Nothing is missing through a lost or corrupted recording.

| File | Column | Rows affected | What `NA` means |
| :--- | :--- | :--- | :--- |
| `steps_index.csv`, `perStep_long_multi.csv` | `gait`, `gait4`, `steadiness` | 1224 of 3842 steps | Not applicable. The step did not pass the quality-control gate, so it was never classified. `qcPass` is 0 and `exclusionReason` reads `did not pass the QC gate` on every one of these rows. |
| `perStep_long_multi.csv` | `contactTime`, `contactTime_n` | 76 of 3842 steps | Not available. No foot-marker contact had its midpoint within the step, so marker-based stance duration could not be measured for it. |
| `perStep_long_multi.csv` | `legCompress_n`, `kLeg_n` | 469 of 3842 steps | Not applicable. The spring-mass inversion returned a compression outside 0 to 0.40 L0 (393 steps), or the step has no marker contact time, an input to the inversion (76 steps), so no leg stiffness is defined. |
| `perStep_long_multi.csv`, `perStride_long_multi.csv` | `dutyLimb` | 964 of 3842 steps, 329 of 1678 strides | Not available. No single-plate contact (20 ms to 0.6 s, above the larger of 10% of body weight and 20% of the plate's peak force) has its midpoint within the step, so the plate-based duty-factor cross-check has no contact to measure. Most fall in the RVC five-plate trials. |
| `perStep_long_multi.csv` | `gaitLabelHand` | all 3842 steps | Not available. The hand gait labels of the original single-study collections were not carried into the compiled roster this package is built from. The column is kept because the analysis code reads it, and the analysis does not use its values. |
| `perStride_long_multi.csv` | `gaitLabelHand`, `gaitCodeHand` | all 1678 strides | As above. |
| `trialRoster_multi.csv` | `logFile` | 266 of 334 trials | Not applicable. The trial has no acquisition log file; only the 68 `Blum_Surface_2012-06` trials were recorded with one. The plate count recorded in a log file is in `trial_index.csv` for every trial, and the sampling rates are in each per-trial `.mat` (`fs_force_Hz`, `fs_kin_Hz`). |

Every other column has a value in every row. `dutyFactor`, the marker-based per-limb duty factor
the paper reports, is defined for all 3842 steps, as is the dimensionless speed `u`. The 334
per-trial time-series CSVs and `mean_cycles.csv` contain no missing values at all. `fracEG` and `CoTmech` divide by step (stride) length and by mean speed, and hold placeholder values of magnitude 10¹² or more, of either sign in `fracEG`, in the 18 steps and 5 strides where that length or speed is not positive; none of these rows is in the analysis sample.

The `.mat` files hold the same content in MATLAB types. In the `steps` table
inside each per-trial `.mat`, the three label columns are empty strings ("") on the 1224 rows
where the CSV reads `NA`, which is how MATLAB represents an unset string, and numeric quantities
use `NaN`.

---

## per_trial_timeseries Folder

Whole on-plate bout per trial on the 500 Hz force time base, as a CSV and a MATLAB `.mat` pair.
All values are in SI units and in the bird-centered frame described under Coordinate Frame below.

### Trajectory File Columns (`.csv`)

| Column        | Units | Description                                                             |
| :------------ | :---- | :----------------------------------------------------------------------- |
| `t_s`         | s     | Time from the start of the analysis window                              |
| `F_ml_N`      | N     | Medio-lateral net ground reaction force, 50 Hz low-pass filtered        |
| `F_fa_N`      | N     | Fore-aft net ground reaction force, positive in the direction of travel |
| `F_vt_N`      | N     | Vertical net ground reaction force, positive up                         |
| `com_ml_m`    | m     | Reconstructed CoM medio-lateral position                                |
| `com_fa_m`    | m     | Reconstructed CoM fore-aft position                                     |
| `com_vt_m`    | m     | Reconstructed CoM vertical position                                     |
| `comProxy_ml_m` | m   | Kinematic CoM proxy, medio-lateral: the cranial and caudal back-marker midpoint with the per-session CoM offset applied (see Reconstruction method) |
| `comProxy_fa_m` | m   | Kinematic CoM proxy, fore-aft                                           |
| `comProxy_vt_m` | m   | Kinematic CoM proxy, vertical                                           |
| `comVel_ml_ms`  | m/s | CoM medio-lateral velocity                                              |
| `comVel_fa_ms`  | m/s | CoM fore-aft velocity                                                   |
| `comVel_vt_ms`  | m/s | CoM vertical velocity                                                   |
| `footR_ml_m`, `footR_fa_m`, `footR_vt_m` | m | Right foot: midpoint of the tarsometatarsophalangeal and third-digit markers, with spikes removed in all three axes together, gaps filled and a 20 Hz low-pass filter |
| `footL_ml_m`, `footL_fa_m`, `footL_vt_m` | m | Left foot, as for the right                                    |
| `stepIndex`   | -     | Step number each sample belongs to, 0 before the first step and after the last; a step's `endSample` is also the next step's `startSample` and is numbered with the next step; joins to `steps_index.csv` |

The reconstructed CoM and the proxy are both provided so a reader can see the reconstruction quality directly for any trial.

The foot side labels come from the marker names in the capture files and were checked against position, since the right foot lies at smaller medio-lateral position than the left (see Coordinate Frame). In 41 trials the two labels were exchanged for the whole capture, and this version exchanges the tracks back; in 24 trials the side cannot be confirmed for the whole trial, and in 4 trials a track does not follow a foot. `foot_side_corrections.csv` lists all 69. The analysis does not depend on the labels, because it takes the stance foot as the lower foot.

The zero-phase 50 Hz filter rings on both sides of a sharp change in force, so `F_vt_N` dips below zero in aerial phases: over the aerial-run steps of the analysis sample, the median of each step's minimum is -0.05 body weights.

### MAT File Contents

Each `.mat` holds one `trial` struct with the same channels as the CSV, as arrays named `time_s` and `stepIndexPerSample` (one column each) and `GRF_N`, `com_m`, `comProxy_m`, `comVel_ms`, `footR_m` and `footL_m` (three columns each, in [ml, fa, vt] order), plus:

* **Acquisition metadata**: `boutID`, `bird`, `subjectID`, `study`, `dateCode`, `colourCode`, `mass_kg`, `L0_m`, `bodyWeight_N`, `fs_force_Hz`, `fs_kin_Hz`, `g`, `nPlates`, `foreaftAxis`, `axisOrder`
* **Reconstruction quality**: `driftRMS_vert_mm`, `WE_identity_r`
* **`steps`**: a table of this trial's rows of `steps_index.csv`, with the same 14 columns

---

## trial_index.csv

One row per trial (334 rows, 15 columns).

| Column             | Units | Description                                                                       |
| :----------------- | :---- | :--------------------------------------------------------------------------------- |
| `boutID`           | -     | Trial identifier, the join key across the package                                 |
| `study`            | -     | Study tag (see Studies Included)                                                  |
| `dateCode`         | -     | Recording date, YYYYMMDD                                                          |
| `bird`             | -     | Study-session bird code                                                           |
| `colourCode`       | -     | Colour tag within the collection cohort                                           |
| `subjectID`        | -     | Reconciled individual identity across cohorts; use this for grouping              |
| `mass_kg`          | kg    | Body mass for that session, estimated from force (see Subject Information)        |
| `L0_m`             | m     | Leg length for that session                                                       |
| `nPlates`          | -     | Number of active force plates                                                     |
| `nSamples`         | -     | Samples in the analysis window                                                    |
| `duration_s`       | s     | Duration of the analysis window                                                   |
| `nSteps`           | -     | Steps detected in the trial                                                       |
| `driftRMS_vert_mm` | mm    | RMS vertical drift of the reconstructed CoM against the kinematic proxy           |
| `WE_identity_r`    | -     | Pearson correlation between the cumulative work of the ground reaction force and the change in CoM mechanical energy over the bout |
| `foreaftAxis`      | -     | Which raw plate axis was resolved as anatomical fore-aft                          |

---

## steps_index.csv

One row per step (3842 rows, 14 columns). The compact step-level index; the full measures are in `perStep_long_multi.csv`.

| Column            | Units | Description                                                                          |
| :---------------- | :---- | :------------------------------------------------------------------------------------ |
| `boutID`          | -     | Trial identifier                                                                     |
| `stepIndex`       | -     | Step number within the trial                                                         |
| `startSample`     | -     | First sample of the step, indexing the per-trial time series                         |
| `endSample`       | -     | Last sample of the step, which is also the next step's `startSample`                 |
| `startTime_s`     | s     | Time of step start                                                                   |
| `endTime_s`       | s     | Time of step end                                                                     |
| `gait`            | -     | `walk`, `groundedRun` or `aerialRun`; the classification the paper uses              |
| `gait4`           | -     | `walkGrounded`, `groundedRun`, `aerialRun` or `pendularRun`; as `gait`, with the aerial-pendular cell (the pendular run) split out of aerial running |
| `steadiness`      | -     | `steady` where the energy grade \|ΔE_CoM/(mgL)\| ≤ 0.05, otherwise `accelerating` or `decelerating` by its sign. Graded on the stride (L the stride length) and inherited by both its steps where the stride passed the stride gate, otherwise graded on the step itself (L the step length; 536 of the 2578 analysis-sample steps) |
| `outlier`         | -     | Mis-cut cycle flag                                                                   |
| `qcPass`          | -     | Retention-gate flag                                                                  |
| `analysisSample`  | -     | 1 for the rows behind the paper                                                      |
| `exclusionReason` | -     | Why a row is not in the analysis sample                                              |
| `u`               | -     | Dimensionless speed, v/√(gL₀)                                                 |

---

## perStep_long_multi.csv

The tidy per-step analysis table: 3842 rows, 90 columns, one row per step. This and `perStride_long_multi.csv` are the tables the analysis and the summary statistics are computed from.

Columns ending `_n` are dimensionless (see Normalisation). Raw and dimensionless forms are both kept so a reader can work in either.

**Identity and grouping**

| Column | Description |
| :--- | :--- |
| `boutID`, `stepIndex` | Trial and step, the join key to `steps_index.csv` and the per-trial series |
| `session`, `bird`, `colourCode`, `study`, `dataset`, `dateCode`, `trialNo` | Session, cohort and collection labels; `session` is `bird\|YYYY-MM`, the grain at which mass and leg length are set |
| `mass_kg`, `L0_m`, `BW_N` | Per-session body mass, leg length and body weight, the normalisation scales |
| `gaitLabelHand` | Legacy hand label from the original collection; not used in the analysis |

**Classification and sample flags**

| Column | Description |
| :--- | :--- |
| `gait`, `gait4` | The classification the paper uses: aerial phase (`hasFlight`) crossed with the direction of energy exchange, pendular where `congruity` < 50% |
| `gaitHodo`, `mechHodo` | The alternative classification by the sign of the CoM velocity-loop area, which the paper tests against `gait`: `mechHodo` is the sign alone (A ≥ 0 pendular), and `gaitHodo` is `aerialRun` wherever `hasFlight` = 1 and otherwise walk or grounded run by that sign |
| `steadiness`, `accClass`, `accClass_fa` | Steadiness class by the energy grade (analysis criterion), by each row's own energy grade, and by the fore-aft speed change |
| `analysisSample`, `exclusionReason`, `outlier` | Sample flags, as in `steps_index.csv` |

**Timing and duty**

| Column | Units | Description |
| :--- | :--- | :--- |
| `stepPeriod`, `stepPeriod_n` | s | Step duration |
| `stepFreq`, `stepFreq_n` | Hz | Step frequency |
| `contactTime`, `contactTime_n` | s | Marker-based stance duration: the longest foot-marker contact of either foot whose midpoint falls within the step (foot height below its 5th percentile plus 3 cm, for at least 40 ms) |
| `flightTime`, `flightTime_n` | s | Time within the step with summed vertical force below 10% of body weight (30 Hz low-pass) |
| `hasFlight` | - | 1 where the summed vertical force is below 10% of body weight for more than 10 ms in total within the step (`flightTime` above 0.010 s) |
| `dutyFactor` | - | Per-limb duty factor from foot-marker contact: the fraction of samples with the foot in contact, averaged over both feet, over the two-step window starting at the step (the last step of a trial uses itself and the step before); the duty factor the paper reports, defined for every step |
| `dutyFactor_gf` | - | Force-based estimate, grounded fraction divided by two, kept for cross-check |
| `groundedFrac` | - | Fraction of the step with summed vertical force above 10% of body weight |
| `dutyLimb` | - | Plate-based duty-factor cross-check: the longest single-plate contact whose midpoint lies within the step, divided by the stride period; above 1 in 27 steps where one plate records a contact longer than the stride period. Not used in the analysis |

**Speed and distance**

| Column | Units | Description |
| :--- | :--- | :--- |
| `meanSpeed`, `meanSpeed_n` | m/s | Mean forward speed over the step. `meanSpeed_n` is the dimensionless speed u |
| `speed_TD`, `speed_TD_n` | m/s | Forward speed at touchdown |
| `stepLength`, `stepLength_n` | m | Fore-aft CoM displacement over the step |
| `Froude` | - | Froude number, u² |
| `accelForeAft`, `accelForeAft_n`, `accelSign` | m/s² | Mean fore-aft acceleration and its sign |
| `vertExcursion`, `vertExcursion_n` | m | Vertical CoM excursion |

**Forces and impulse**

| Column | Units | Description |
| :--- | :--- | :--- |
| `peakVertForce_BW` | BW | Peak vertical ground reaction force |
| `peakForeAftForce_BW` | BW | Peak absolute fore-aft force |
| `peakResultantForce_BW` | BW | Peak resultant force magnitude |
| `impulseForeAft`, `impulseForeAft_n` | N s | Net fore-aft impulse over the step, touchdown to next touchdown |

**Center-of-mass energy and work**

| Column | Units | Description |
| :--- | :--- | :--- |
| `dKE_step`, `dKE_step_n` | J | Change in total CoM kinetic energy ½m(v_ml² + v_fa² + v_vert²) over the step |
| `dPE_step`, `dPE_step_n` | J | Change in gravitational potential energy mgz over the step; `dKE_step` + `dPE_step` = `dE_CoM` |
| `dE_CoM`, `dE_CoM_n` | J | Net change in total CoM mechanical energy |
| `W_CoM_pos`, `W_CoM_neg`, `W_CoM_net` (and `_n`) | J | Positive, negative and net CoM work |
| `KE_range`, `PE_range` (and `_n`) | J | Within-step range of the forward kinetic energy E_kf = ½m(v_ml² + v_fa²) and of mgz |
| `ampRatio_KEtoPE` | - | Range of E_kf divided by the range of the vertical energy E_v = mgz + ½mv_vert² |
| `WE_relresid`, `WE_r` | - | `WE_relresid`: \|W_CoM_net - ΔE_CoM\| divided by the within-step range of CoM mechanical energy (in the stride table, by W_f + W_v); the quality-control gate requires at most 0.20. `WE_r`: the trial's work-energy correlation `WE_identity_r`, repeated on each row |
| `fracEG`, `fracDE`, `fracDV` | - | Steadiness measures: energy grade ΔE_CoM/(mgL), with L the step length (the stride length in the stride table), the analysis criterion; energy fraction ΔE_CoM/(mv²); and the signed fore-aft speed change aT/v |
| `driftRMS_vert_mm` | mm | Trial-level vertical reconstruction drift, repeated on each step |

**Gait descriptors**. The analysis sample requires six descriptors to be defined: `recovery`, `congruity`, `CoTmech`, `collisionAngle`, `dutyFactor` and the touchdown mechanical energy, which `03_step_descriptors.R` computes from `legLen_TD_n`, `legAngle_TD` and the touchdown velocity and which is not deposited. `collisionFraction` and the loop area are given as well.

| Column | Units | Description |
| :--- | :--- | :--- |
| `recovery` | % | Pendular energy recovery over the step (Cavagna, Heglund and Taylor 1977) |
| `congruity` | % | Percentage of the step over which E_kf and E_v change in the same direction (Ahn, Furrow and Biewener 2004) |
| `collisionAngle` | rad | Force- and velocity-weighted deviation of the ground reaction force from perpendicular to CoM velocity (Lee et al. 2011, 2013) |
| `collisionFraction` | - | Force- and velocity-weighted mean of φ/(θ + λ), the collision angle relative to its potential, with θ the deviation of the force from vertical and λ that of the velocity from horizontal (Lee et al. 2011, 2013) |
| `CoTmech` | - | Dimensionless mechanical cost of transport |
| `hodoArea`, `hodoArea_n` | m²/s² | Signed area of the CoM velocity loop by the shoelace formula. Positive is counterclockwise (pendular), negative clockwise (bouncing). `hodoArea_n` is A/(gL₀) |

**Virtual leg and leg stiffness**

| Column | Units | Description |
| :--- | :--- | :--- |
| `legLen_TD`, `legLen_TD_n` | m | CoM-to-foot distance at touchdown |
| `legAngle_TD` | deg | Virtual leg angle at touchdown, counterclockwise from horizontal with fore-aft positive in the direction of travel |
| `dLegLen`, `dLegAngle` | m, deg | Change in leg length and angle over stance; a negative `dLegAngle` is limb retraction |
| `legCompress_n` | - | Leg compression as a fraction of the session leg length `L0_m`, which is also the leg length the inversion uses |
| `kLeg_n` | BW/L0 | Dimensionless leg-spring stiffness kL₀/(mg), from k = F_max/ΔL with compression recovered by inverting a sinusoidal vertical force profile over contact (Blum, Lipfert and Seyfarth 2009, method C) |

`legCompress_n` and `kLeg_n` are `NA` where the inversion returns a compression outside 0 to 0.40 L0, or where `contactTime` is `NA`. Leg stiffness is defined only relative to the spring-mass template it comes from, and where two limbs share the net force in double support it is an effective two-limb stiffness rather than a single-limb property.

---

## perStride_long_multi.csv

The tidy per-stride table: 1678 rows, 78 columns, one row per stride. The stride with `strideIndex` k spans steps k and k + 1 of the same trial, for odd k (see Data Selection). Column meanings match `perStep_long_multi.csv` except as noted.

| Column group | Columns | Notes |
| :--- | :--- | :--- |
| Identity | `boutID`, `strideIndex`, `session`, `bird`, `colourCode`, `study`, `dataset`, `dateCode`, `trialNo`, `mass_kg`, `L0_m`, `BW_N`, `gaitLabelHand`, `gaitCodeHand` | `strideIndex` is the `stepIndex` of the stride's first step; its second step is `strideIndex` + 1 |
| Sample flag | `analysisSample`, `outlier` | The stride table carries no `gait` column; join to the per-step table for gait |
| Timing | `stridePeriod`, `strideFreq`, `contactTime`, `flightTime`, `hasFlight`, `dutyFactor`, `dutyFactor_gf`, `groundedFrac`, `dutyLimb` (and `_n` forms) | Means of the stride's two step values; `stridePeriod` and `strideFreq` are taken over the stride, and `hasFlight` is 1 where either step has an aerial phase |
| Speed and distance | `meanSpeed`, `strideLength`, `Froude`, `accelForeAft`, `accelSign`, `vertExcursion` (and `_n` forms) | |
| Energy and work | `dKE`, `dPE`, `Wf_horiz`, `Wv_vert`, `Wtot_com`, `W_CoM_pos`, `W_CoM_neg`, `W_CoM_net`, `dE_CoM`, `KE_range`, `PE_range`, `ampRatio_KEtoPE`, `WE_relresid`, `WE_r`, `driftRMS_vert_mm` (and `_n` forms) | `Wf_horiz`, `Wv_vert` and `Wtot_com` are the summed positive increments of the horizontal, vertical and total energy that enter the recovery calculation; `dKE` is the change in E_kf and `dPE` the change in mgz over the stride, so they do not sum to `dE_CoM`, which has no `_n` form here |
| Steadiness | `fracEG`, `fracDE`, `fracDV`, `accClass`, `accClass_fa` | Graded on the stride; the per-step `steadiness` inherits it where the stride passed the stride gate |
| Descriptors | `recovery`, `congruity`, `collisionAngle`, `collisionFraction`, `CoTmech`, `hodoArea`, `hodoArea_n`, `gaitHodo`, `mechHodo` | |

Pendular recovery is R = (W_f + W_v - W_tot)/(W_f + W_v), each W the sum of positive increments over the row's own window. Recovery and congruity are computed per step in `perStep_long_multi.csv`, the values the paper uses, and per stride here.

---

## morphology_multi.csv

One row per study session (30 rows, 10 columns).

| Column     | Units | Description                                                              |
| :--------- | :---- | :------------------------------------------------------------------------ |
| `session`  | -     | Session key, `bird\|YYYY-MM`; the grain at which mass and leg length are set |
| `bird`     | -     | Study-session bird code                                                  |
| `study`    | -     | Study tag                                                                |
| `sessionMonth` | - | Month of the recording session, `YYYY-MM`                                 |
| `mass_kg`  | kg    | Body mass for that session, estimated from force (see Subject Information) |
| `L0_m`     | m     | Leg length for that session; `L0_source` gives how it was set           |
| `Lnorm_m`  | m     | Geometric length scale, 0.20 m^(1/3), reported for reference; it is not the normaliser used |
| `nSteps`   | -     | Steps detected for that session                                          |
| `nStrides` | -     | Strides detected for that session                                        |
| `L0_source`| -     | How `L0_m` was set for that session: `session_nonAerial_median`, `bird_nonAerial_median` or `withinSample_allometry` |

---

## trialRoster_multi.csv

Source-file provenance, one row per trial (334 rows, 12 columns): `boutID`, `bird`, `colourCode`, `study`, `dateCode`, `trialNo`, `trialStem`, `bodyMass_kg`, `windowMode`, and the `forceFile`, `markerFile` and `logFile` paths of the raw recordings.

Paths are relative to the root of the source archive they were collated in, for example
`OtherSpecies_byDate/Guinea fowl/Blum_Surface_2012-06/2012-06-25/gf_red_LE0cm_0006.txt`, so they carry no machine or user name. `logFile` is `NA` for the trials that have none: only the 68 `Blum_Surface_2012-06` trials were recorded with one.

The raw force and marker files are not part of this deposit and are not needed to reproduce the analysis, which runs from `per_trial_timeseries/`. This table is provenance: it records which source recording each trial came from, so a trial here can be traced back to the collection it belongs to.

---

## force_sign_corrections.csv

One row per corrected force column (boutID, axis, column, reason): the trials and the axis, `ml` or `fa`, whose ground reaction force had the opposite sign in the first version. Version History describes the correction.

---

## foot_side_corrections.csv

One row per trial the foot-marker side check flags (69 rows): `boutID`, `study`, `correction`, `detail` and `analysisSampleSteps`. `correction` is `swapped` (41 trials: the labels were exchanged for the whole capture, and this version exchanges the tracks back), `inconsistent` (24 trials: the side cannot be confirmed for the whole trial, and the labels are as captured) or `invalidTrack` (4 trials: a foot track's median horizontal distance from the CoM proxy exceeds 0.20 m, so it does not follow a foot; none of their steps is in the analysis sample). `detail` gives the evidence, computed from the tracks as captured, before any exchange: the median of `footR_ml_m` - `footL_ml_m` in each third of the trial, or the track's distance from the CoM proxy.

---

## mean_stepcycle_by_gait_speed

Ensemble mean step cycles for each gait by speed-bin group, in long format (10800 rows) and as a MATLAB struct. Six groups: walk, grounded run and aerial run, each split slow and fast.

Two things determine which cycles enter a group, which corresponds to the paper's mean-trace figure. 1) The cycles are the STEADY steps of the analysis sample (`steadiness` = `steady` and `analysisSample` = 1). 2) The slow/fast boundary is the gait's median Froude number over its analysis-sample steps of every steadiness class, not the median of the steady subset, so the boundary does not move between steadiness classes; the two halves of a group can therefore be unequal in size.

| Column      | Units | Description                                                              |
| :---------- | :---- | :------------------------------------------------------------------------ |
| `gait`      | -     | `walk`, `groundedRun` or `aerialRun`                                     |
| `speedbin`  | -     | `slow` or `fast`, split at the gait's median Froude over its analysis-sample steps |
| `n`         | -     | Number of step cycles in the group mean                                  |
| `dur_mean_s`| s     | Mean cycle duration for the group                                        |
| `pct`       | %     | Percent of the step cycle, 0 to 100                                      |
| `time_s`    | s     | Real-time axis, `pct` scaled by `dur_mean_s`, preserving the mean duration |
| `channel`   | -     | Which of the 18 channels this row belongs to                             |
| `mean`      | varies| Group mean of that channel at that phase                                 |
| `sd`        | varies| Standard deviation across cycles                                         |

Channels: `F_vt_N`, `F_fa_N`, `F_ml_N`, `F_vt_BW`, `F_fa_BW` (forces, raw and body-weight normalised); `vfa_ms`, `vvert_ms` (CoM velocity); `com_fa_disp_m`, `com_vt_rel_m` (CoM displacement, fore-aft from the step start and vertical about the step mean); `KE_J`, `PE_J`, `Etot_J`, `KE_n`, `PE_n` (fluctuations of total kinetic, potential and total CoM energy about each step's mean, raw and normalised); `stanceFoot_fa_relCoM_m`, `stanceFoot_vt_relCoM_m`, `swingFoot_fa_relCoM_m`, `swingFoot_vt_relCoM_m` (foot position relative to the CoM).

---

## Coordinate Frame and Units

Every three-vector is ordered `[medio-lateral, fore-aft, vertical]`. The frame is bird-centered: fore-aft is positive in the direction of travel for every trial, vertical is positive up, and medio-lateral is positive to the bird's left, so fore-aft, medio-lateral and vertical, in that order, form a right-handed frame. Trials in which the bird crossed the runway in the opposite lab direction have their fore-aft and medio-lateral signs flipped together, so all trials are directly comparable without further reorientation.

Units are SI: force in newtons, position in meters, velocity in meters per second, time in seconds, mass in kilograms, energy in joules. Angles are in degrees except the collision angle, which is in radians because the equivalence with the dimensionless cost of transport holds in that form.

## Normalisation

Data are normalised into dimensionless quantities following dynamic similarity, using each session's body mass m, leg length L₀ and g = 9.81 m/s². Columns ending `_n` carry the dimensionless form.

| Quantity | Normalisation factor | Formula        |
| :------- | :------------------- | :------------- |
| Force    | Body weight          | mg           |
| Length   | Leg length           | L₀          |
| Velocity | Froude velocity      | √(gL₀)  |
| Time     | Pendular period      | √(L₀/g) |
| Frequency| Inverse pendular period | √(g/L₀) |
| Work, energy | Body weight times leg length | mgL₀ |
| Impulse  | Body weight times pendular period | mg√(L₀/g) |
| Loop area| Gravity times leg length | gL₀      |
| Acceleration | Gravity | g |

Dimensionless speed is u = v/√(gL₀) and the Froude number is u².

---

## Processing Workflow

The analysis code is on GitHub at [https://github.com/NeuroMechLab/GuineaFowlGait-code](https://github.com/NeuroMechLab/GuineaFowlGait-code), tagged `v1.3.1`. `PIPELINE.md` there gives the full run order and maps every output to the manuscript element that uses it.

### Reproducing the analysis from this package

The release runs from this package. The per-trial files carry the net ground reaction force, the time base, the kinematic CoM proxy the reconstruction is path-matched to, both foot marker tracks, and body mass, leg length and g, which is everything the measurement kernels read, so the raw force-plate and marker recordings are not needed.

Unzip `per_trial_timeseries.zip` and `mean_stepcycle_by_gait_speed.zip` inside the package, put the package at `GaitSel_DryadPackage_AllGF/` beside `_MatlabProcessingCode/`, then:

    matlab -batch "addpath(genpath('_MatlabProcessingCode')); \
        GaitSelMulti_BatchFromDryad; GaitSelMulti_ExportTidyCSV; \
        GaitSelMulti_QCGate; GaitSelMulti_ExportCycleTraces"
    cd _RAnalysis && Rscript run_all.R

**Step 1.** `GaitSelMulti_BatchFromDryad.m` rebuilds each trial from `per_trial_timeseries/` and re-runs the reconstruction (`GaitSel_ReconstructCoM`), the gait-event detection at foot-marker touchdowns (`GaitSel_DetectGaitEvents` with `GaitSel_FootContacts`) and the per-step and per-stride reduction (`GaitSel_PerStepStrideMeasures`).

**Step 2.** `GaitSelMulti_ExportTidyCSV.m` derives per-session `L0`, adds the dimensionless and leg-stiffness columns and flags outliers, writing `perStep_long.csv`, `perStride_long.csv` and `morphology.csv`.

**Step 3.** `GaitSelMulti_QCGate.m` applies the work-energy and drift criteria and the energy-grade steadiness criterion. This is the single definition of the retained set, computed once in MATLAB so that this package and the statistics agree.

**Step 4.** `GaitSelMulti_ExportCycleTraces.m` writes the per-step cycle traces the figures read, also from `per_trial_timeseries/`.

**Step 5.** `Rscript run_all.R` runs the numbered R chain: `01_load.R` reads the tidy tables and builds the individual grouping, `02_clean.R` applies the quality-control and steadiness flags, `03_step_descriptors.R` adds the collision angle and touchdown energy, `04_analysis_sample.R` fixes the analysis sample of 2578 steps and 942 strides, and the remaining scripts produce every figure, table and statistic. `GaitSel_ExportDryadMeanCycles.m` rebuilds the mean step cycles afterwards, since it reads the R phase's speed-bin thresholds, and writes them to `rerun/mean_stepcycle_by_gait_speed/`.

Nothing is written back into this package: the re-run writes to `rerun/`, `SavedResults/` and `_RAnalysis/` (`data/`, `output/`, and `output_internal/`, whose `dryad_tables/` receives the two relabelled tidy tables), so a second run starts from the archived state. `rerun/reconstruction_check.csv` compares the re-run against the values these files already carry. Over all 334 trials it agrees to 0.000 mm on the vertical CoM path and on the drift, to within 3e-8 on the work-energy correlation, and on the number of steps cut for every trial.

Two columns cannot be reproduced from this package. `dutyLimb`, the plate-based duty-factor cross-check, is estimated from the per-plate vertical force, which is not deposited, and `trialNo`, the last field of `boutID`, is not stored in the per-trial series; both come back as `NA`. The duty factor the paper reports is the marker-based `dutyFactor` and is unaffected.

### How this package was made

The per-trial files themselves come from the raw force-plate and marker recordings, which are not part of this deposit. Those were paired into a trial roster, imported and windowed, the horizontal force axes resolved to anatomical fore-aft and medio-lateral, the CoM offset applied to the marker proxy, the foot tracks de-spiked and their side labels checked against position, before the reconstruction and measurement steps above. That import layer has no input a reader can supply and is not part of the code release.

### Reconstruction method

The CoM is reconstructed by path-matched double integration of the ground reaction force over the whole on-plate bout: F = ma, with one set of integration constants per axis (initial position, initial velocity and an acceleration-baseline offset) fitted to minimise departure from the kinematic CoM proxy. The net force is low-pass filtered at 50 Hz with a zero-phase Butterworth filter. The analysis window for each trial is found automatically as the span where the bird is both tracked and on the plates. The orientation of the force plates relative to the motion-capture frame is fixed within a collection date and differs between collection years, so the horizontal force axes are assigned once per date. Fore-aft is the plate axis, with the sign, whose cumulative impulse best follows the kinematic fore-aft velocity over all of that date's trials. Medio-lateral is the other plate axis, with the sign that keeps the plate frame and the motion-capture frame both right-handed.

The proxy is the midpoint of the cranial and caudal back markers, which lies above and behind the whole-body CoM, so a fore-aft and a vertical offset are applied to it before any geometry is taken. The fore-aft offset is fitted by driving the net pitch impulse about the CoM to zero over stride-complete windows, within ±0.60 L_iso of the back-marker midpoint; no session's fit reached that bound. That condition does not identify the vertical offset, because the fore-aft impulse of the summed force approaches zero over a whole stride, so the vertical offset is set from anatomical scaling as -0.2625 L_iso with L_iso = 0.20 m^(1/3). Offsets are fitted for each bird in each recording session, because marker placement was renewed at each session. `comProxy_*` in the per-trial series carries them, so leg length, leg angle and CoM height above the foot all refer to one CoM. Energy fluctuations are unaffected, because a constant shift cancels in a height change.

Reconstruction quality per study, as median work-energy correlation r and median vertical drift: RVC 0.99 and 7.7 mm, DropVsPothole 1.00 and 7.5 mm, Surface 0.97 and 4.9 mm.

---

## Version History

### Version 3: foot-marker side labels and README corrections

The side labels of the foot tracks (`footR_*`, `footL_*`; `footR_m`, `footL_m` in the .mat files) were checked against position. In 41 trials, 37 of the 77 `Blum_DropVsPothole_2012-02` trials, 3 `RVC_Daley_2008-2011` trials and 1 `Blum_Surface_2012-06` trial, the right and left labels were exchanged for the whole capture, and this version exchanges the two tracks back. The 24 trials whose side cannot be confirmed for the whole trial and the 4 trials with a track that does not follow a foot are listed, unchanged, in `foot_side_corrections.csv`. The analysis takes the stance foot as the lower foot and treats the two feet alike, so no step, stride, label or measure changes: the code release run on this version and on version 2 writes identical output files, and each of the 89 values the article reports is the same as under version 2.

This README also corrects statements that did not match the code or the files:

* the coordinate frame: medio-lateral is positive to the bird's left;
* column definitions: `mass_kg`, `L0_m`, `WE_identity_r`, `WE_r`, `WE_relresid`, `dKE_step`, `dPE_step`, `KE_range`, `PE_range`, `ampRatio_KEtoPE`, `fracEG`, `fracDV`, `gait`, `gaitHodo`, `mechHodo`, `recovery`, `congruity`, `hasFlight`, `flightTime`, `groundedFrac`, `dutyFactor`, `dutyLimb`, `contactTime`, `impulseForeAft`, `collisionFraction`, `legCompress_n`, `steadiness`, `stepIndex`, `endSample` and `strideIndex`, the foot tracks, the stride timing and energy columns, the mean-cycle speed split and energy channels, and the array names in the .mat files;
* selection rules and methods: the leg-length window behind L0, the drift and touchdown-height bounds, the descriptors the analysis sample requires, the causes of missing leg stiffness and `dutyLimb`, where the sampling rates are stored, the placeholder values in `fracEG` and `CoTmech`, the CoM offset bound, the toolboxes used, the trial-name formats and the join keys;
* the re-run: the files it writes, the two columns it cannot reproduce, and the analysis sample it gives (2578 steps and 942 strides);
* the number of empty labels in the .mat step tables (1224).

The pipeline forms no stride from the final pair of steps in a trial with an even number of steps (see Data Selection); this version documents that, and the stride and steadiness columns are as in version 2.

### Version 2: horizontal force sign correction

Version 2 corrected two errors in the sign of the horizontal ground reaction force (`F_ml_N`, `F_fa_N`; `GRF_N` in the .mat files).

* Medio-lateral force, in all 74 trials recorded on 2008-11-20, 2008-11-21, 2009-11-17, 2009-11-18, 2009-11-19. In these sessions the force plates were turned 180 degrees relative to the motion-capture frame, and the medio-lateral axis had been given the orientation of the other sessions.
* Fore-aft force, in 4 trials whose own sign fit disagreed with the other trials of their collection date: `rvc_red_20081120_L00cm_0006`, `surf_noc_20120626_LE0cm_0003`, `surf_red_20120625_LE0cm_0004`, `surf_yel_20120627_LE0cm_0009`. `rvc_red_20081120_L00cm_0006` is one of the four `rvc_red` trials of 2008-11-20 whose force and marker files come from different recordings; its sign follows the other trials of that date, and none of those four trials has a step in the analysis sample.

Integrated over body mass, the listed force components changed the CoM velocity in the direction opposite to the velocity of the CoM proxy, the back-marker midpoint. With the corrected signs, the reconstructed CoM follows the proxy: over the 74 trials with corrected medio-lateral force, the median root-mean-square distance between them in the medio-lateral direction is 7.5 mm, against 27.6 mm in the first version.

Recomputed from the corrected force: the CoM path and velocity, the drift and the work-energy identity of the corrected trials; the leg length L0 of the 12 sessions in which it changes by more than one micrometre, since L0 is a median CoM-to-foot distance (the largest changes are `rvc_yel|2008-11`, 0.195 to 0.163 m, and `rvc_blu|2008-11`, 0.203 to 0.187 m); every quantity normalised by L0 in those sessions; the per-step and per-stride measures, labels and flags of the 95 affected trials; and the mean step cycles. Of the affected trials, 16 lie in the two allometry sessions, whose L0 depends on the other sessions, and 2 change only in stride flags and, for one step, the steadiness label inherited from a re-flagged stride, because the mis-cut-cycle filter compares each stride with strides of similar dimensionless speed in all trials. Every other row and per-trial file is identical to the first version.

In the November 2008 sessions of `rvc_yel` and `rvc_blu`, the corrected leg lengths are shorter than in the same birds' later sessions (`rvc_yel` 0.163 m against 0.197 m in April 2009, `rvc_blu` 0.187 m against 0.200 m in May 2009), because the height of the CoM above the stance foot at touchdown is lower in that collection (0.134 m against 0.180 m for `rvc_yel`), which is likely due to differences in back-marker placement between recording sessions; these two sessions provide 30 of the 2578 analysis-sample steps, so the difference is unlikely to shift any major result of the article.

With version 2 the analysis sample is 2578 steps and 942 strides, against 2580 and 944. Of the 89 values the article reports, 51 change. The three crossover speeds are u = 1.057 (Froude 1.116), u = 0.560 (Froude 0.313) and u = 0.576 (Froude 0.332) (article: u = 1.056 (Froude 1.114), u = 0.557 (Froude 0.310) and u = 0.576 (Froude 0.332)); the agreement between the two-cluster split and rotation sense is 87.7% (2260 of 2578 steps) (article: 87.7% (2262 of 2580 steps)); and the leg-length range is 0.163-0.232 m (article: 0.192-0.232 m).

Although the numbers shift slightly, none of the claims in the paper are changed by this update.  The main claims hold at the precision reported in the paper. Crossover speeds: rotation sense u = 0.58 (0.46 to 0.63) and footfall u = 1.06 (1.01 to 1.09) are unchanged when rounded as in the text. The energy-phase crossover stays at u = 0.56. The rotation-sense and energy-phase intervals still overlap, and the footfall transition is still roughly twice as fast. Agreement with the cluster split: 87.7% unchanged. Descriptive ordering across walk, grounded run and aerial run: each median moves by at most 0.4 units, and the order is the same. Leg stiffness: 7.1, 6.2 and 9.7 for those three gait categories is unchanged.

---

## Code Dependencies

**MATLAB** R2019b or later, with the Statistics and Machine Learning and Signal Processing Toolboxes.

**R** 4.6.0, with `dplyr`, `ggplot2`, `patchwork`, `readr`, `tidyr`, `mgcv`, `cluster`, `openxlsx` and `knitr`.

All scripts locate the project root from their own file location, so no path configuration is required.

---

## License

Data in this deposit are released under **CC0 1.0 Universal** public domain dedication, as Dryad requires: [https://creativecommons.org/publicdomain/zero/1.0/](https://creativecommons.org/publicdomain/zero/1.0/)

The analysis code on GitHub is released separately under the MIT License.

---

## Citation

Please cite the associated manuscript when using this data or code:\
Daley, M. A. and Birn-Jeffery, A. (2026). A velocity-loop view of avian gait: hodograph rotation
sense characterizes the walking-to-running continuum in guinea fowl. *Biology Open*,
doi:10.1242/bio.062880

Please also cite this data package:

Daley, M. A. and Birn-Jeffery, A. (2026). Data from: A velocity-loop view of avian gait:
hodograph rotation sense characterizes the walking-to-running continuum in guinea fowl. Dryad,
doi:10.5061/dryad.7m0cfxqc4

Contact: For questions regarding this dataset or code, please contact the corresponding author Monica Daley ([madaley@uci.edu](mailto:madaley@uci.edu)).
