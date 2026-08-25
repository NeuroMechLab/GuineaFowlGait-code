## Guinea fowl center-of-mass dynamics across the walking-to-running speed range

Data and code associated with the paper: **"A velocity-loop view of avian gait: hodograph rotation sense characterizes the walking-to-running continuum in guinea fowl"** Authors: Monica A. Daley and Aleksandra Birn-Jeffery\
*Biology Open* (2026). doi:10.1242/bio.062880

\***Author for correspondence:** [madaley@uci.edu](mailto:madaley@uci.edu)

## Overview

This package contains ground reaction force recordings and reconstructed center-of-mass (CoM) trajectories for guinea fowl (*Numida meleagris*) running over level ground, compiled from three Royal Veterinary College data collections and reprocessed through one common pipeline. The analysis asks how the rotation sense of the CoM velocity loop, the hodograph, relates to the transition from pendular (walking) to bouncing (running) dynamics, and compares it against established center-of-mass criteria across a speed range that includes extensive grounded running.

Each of 334 trials is provided as a whole-recording time series in a standardized, bird-centered coordinate frame, together with tidy per-step and per-stride measure tables, per-step gait classification, and ensemble mean step cycles by gait and speed.

## Directory Structure

```
GaitSel_DryadPackage_AllGF/
├── README.md                                      # This file
├── trial_index.csv                                # One row per trial (334)
├── steps_index.csv                                # One row per step (3842)
├── perStep_long_multi.csv                         # Tidy per-step measures (3842 x 85)
├── perStride_long_multi.csv                       # Tidy per-stride measures (1678 x 77)
├── morphology_multi.csv                           # Per-session body mass and leg length (30)
├── trialRoster_multi.csv                          # Source file provenance per trial (334)
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
| trialStem | L00cm_####     | Terrain code and sequential trial number within that session       |

The same string is the `boutID`, which is the join key across every table in this package.

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

Thirty bird-sessions over thirteen individuals. Colour codes name the same individual within a collection cohort but were reused across cohorts, so `subjectID` is the reconciled identity and is the correct grouping factor for any analysis pooling across studies. The two 2012 Blum studies, four months apart, are one cohort. The June bird recorded as `noc`, carrying no colour band, is the blue bird: the study's weight record lists blue among the seven birds run and the archive holds no `blu` for that study, it is the only unmatched bird, and the original analysis's published data uses blue's weighed mass for 'noc'. Body mass and leg length are per session, since a bird's mass differs between sessions.

| Session | subjectID | Study | Session date | Mass (kg) | L0 (m) | L0 source | Trials | Steps | Strides |
| :------ | :-------- | :---- | :----------- | --------: | -----: | :-------- | -----: | ----: | ------: |
| `pot_bla` | `blum12_bla` | Blum_DropVsPothole_2012-02 | 2012-02 | 1.158 | 0.208 | session | 10 | 169 | 77 |
| `pot_blu` | `blum12_blu` | Blum_DropVsPothole_2012-02 | 2012-02 | 1.087 | 0.205 | session | 10 | 133 | 60 |
| `pot_g00` | `blum12_g00` | Blum_DropVsPothole_2012-02 | 2012-02 | 1.313 | 0.205 | session | 10 | 132 | 58 |
| `pot_gry` | `blum12_gry` | Blum_DropVsPothole_2012-02 | 2012-02 | 1.377 | 0.222 | session | 10 | 140 | 61 |
| `pot_reb` | `blum12_reb` | Blum_DropVsPothole_2012-02 | 2012-02 | 1.177 | 0.219 | session | 10 | 134 | 60 |
| `pot_red` | `blum12_red` | Blum_DropVsPothole_2012-02 | 2012-02 | 1.243 | 0.204 | session | 10 | 181 | 83 |
| `pot_rre` | `blum12_rre` | Blum_DropVsPothole_2012-02 | 2012-02 | 1.253 | 0.208 | allometry | 7 | 52 | 20 |
| `pot_yel` | `blum12_yel` | Blum_DropVsPothole_2012-02 | 2012-02 | 1.172 | 0.198 | session | 10 | 121 | 53 |
| `surf_bla` | `blum12_bla` | Blum_Surface_2012-06 | 2012-06 | 1.453 | 0.222 | session | 10 | 190 | 87 |
| `surf_g00` | `blum12_g00` | Blum_Surface_2012-06 | 2012-06 | 1.740 | 0.219 | session | 9 | 131 | 59 |
| `surf_gry` | `blum12_gry` | Blum_Surface_2012-06 | 2012-06 | 1.492 | 0.232 | session | 10 | 111 | 48 |
| `surf_noc` | `blum12_blu` | Blum_Surface_2012-06 | 2012-06 | 1.361 | 0.210 | session | 10 | 132 | 59 |
| `surf_red` | `blum12_red` | Blum_Surface_2012-06 | 2012-06 | 1.332 | 0.214 | session | 10 | 129 | 58 |
| `surf_rre` | `blum12_rre` | Blum_Surface_2012-06 | 2012-06 | 1.424 | 0.207 | allometry | 9 | 68 | 26 |
| `surf_yel` | `blum12_yel` | Blum_Surface_2012-06 | 2012-06 | 1.302 | 0.220 | session | 10 | 106 | 45 |
| `rvc_blu` | `rvc_blu` | RVC_Daley_2008-2011 | 2008-11 | 1.161 | 0.203 | session | 2 | 18 | 8 |
| `rvc_blu` | `rvc_blu` | RVC_Daley_2008-2011 | 2009-05 | 1.378 | 0.200 | session | 13 | 113 | 48 |
| `rvc_blu` | `rvc_blu` | RVC_Daley_2008-2011 | 2009-11 | 1.286 | 0.203 | bird | 10 | 47 | 16 |
| `rvc_bly` | `rvc_bly` | RVC_Daley_2008-2011 | 2008-11 | 1.277 | 0.197 | bird | 4 | 14 | 4 |
| `rvc_bly` | `rvc_bly` | RVC_Daley_2008-2011 | 2009-05 | 1.653 | 0.197 | session | 29 | 471 | 214 |
| `rvc_bly` | `rvc_bly` | RVC_Daley_2008-2011 | 2009-11 | 1.273 | 0.197 | bird | 11 | 80 | 32 |
| `rvc_red` | `rvc_red` | RVC_Daley_2008-2011 | 2008-11 | 1.380 | 0.216 | session | 4 | 35 | 14 |
| `rvc_red` | `rvc_red` | RVC_Daley_2008-2011 | 2009-05 | 1.934 | 0.197 | session | 20 | 286 | 128 |
| `rvc_red` | `rvc_red` | RVC_Daley_2008-2011 | 2009-11 | 1.338 | 0.200 | bird | 12 | 86 | 34 |
| `rvc_rey` | `rvc_rey` | RVC_Daley_2008-2011 | 2008-11 | 1.317 | 0.193 | session | 1 | 9 | 4 |
| `rvc_rey` | `rvc_rey` | RVC_Daley_2008-2011 | 2009-05 | 1.771 | 0.192 | session | 24 | 363 | 161 |
| `rvc_rey` | `rvc_rey` | RVC_Daley_2008-2011 | 2009-11 | 1.276 | 0.208 | session | 15 | 107 | 44 |
| `rvc_yel` | `rvc_yel` | RVC_Daley_2008-2011 | 2008-11 | 1.203 | 0.195 | session | 4 | 31 | 13 |
| `rvc_yel` | `rvc_yel` | RVC_Daley_2008-2011 | 2009-04 | 1.504 | 0.198 | session | 29 | 313 | 135 |
| `rvc_yel` | `rvc_yel` | RVC_Daley_2008-2011 | 2009-11 | 1.190 | 0.198 | bird | 11 | 59 | 22 |

Leg length `L0` is the median CoM-to-foot distance at touchdown over the session's steps with no aerial phase. Where a session yielded fewer than five such steps it falls back to that bird's sessions pooled, and for the one bird that never yielded five, to a within-sample allometry on body mass; the `L0 source` column records which applies.

### Data Selection

Every detected step and stride is provided, with flags rather than deletion, so a reader can reproduce the paper's sample or define a different one.

```
3842 steps detected
2620 pass the quality-control gate           (qcPass = 1)
  77 duty factor undetermined
  47 touchdown height out of range
   1 both
2580 in the analysis sample                  (analysisSample = 1)
```

| Flag             | Meaning                                                                                                                                                                                                                                     |
| :--------------- | :------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `outlier`        | 1 marks a mis-cut cycle: a log-MAD filter on step frequency and length judged within dimensionless-speed neighborhoods, combined with a speed-conditional filter catching merged (half-frequency) and split (double-frequency) cycles       |
| `qcPass`         | 1 marks steps passing the retention gate: `outlier = 0`, finite positive speed, work-energy identity residual at or below 0.20, and vertical reconstruction drift within a duration-normalized bound (the larger of 45 mm and a fixed per-second rate times trial duration) |
| `analysisSample` | 1 marks the rows behind every number in the paper. Adds to the retention gate the requirement that every continuous gait-space descriptor be defined                                                                                        |
| `exclusionReason`| Names the reason for every row outside the analysis sample                                                                                                                                                                                  |

Two gait-space descriptors depend on foot-marker quality that cannot be detected based on a force-based gate: the per-limb duty factor is undefined where foot contact could not be determined from the markers, and the touchdown energy is considered infeasible where the reconstructed CoM height above the contacting foot falls outside 0.4 to 1.6 leg lengths.

**`analysisSample = 1` is the column to filter on to reproduce the paper.** It selects 2580 steps and 944 strides. Within the step sample the steadiness split is 930 steady, 947 accelerating and 703 decelerating. Of 1678 strides detected, 1129 pass the stride gate and 944 are in the analysis sample, those being the strides that pass the gate and whose two steps are both in the step sample.

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
| `footR_ml_m`, `footR_fa_m`, `footR_vt_m` | m | Right foot marker position, de-spiked and gap-filled     |
| `footL_ml_m`, `footL_fa_m`, `footL_vt_m` | m | Left foot marker position, de-spiked and gap-filled      |
| `stepIndex`   | -     | Step number each sample belongs to; joins to `steps_index.csv`          |

The reconstructed CoM and the proxy are both provided so a reader can see the reconstruction quality directly for any trial.

### MAT File Contents

Each `.mat` holds one `trial` struct with the same channels as the CSV plus:

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
| `mass_kg`          | kg    | Body mass for that session                                                        |
| `L0_m`             | m     | Leg length for that session                                                       |
| `nPlates`          | -     | Number of active force plates                                                     |
| `nSamples`         | -     | Samples in the analysis window                                                    |
| `duration_s`       | s     | Duration of the analysis window                                                   |
| `nSteps`           | -     | Steps detected in the trial                                                       |
| `driftRMS_vert_mm` | mm    | RMS vertical drift of the reconstructed CoM against the kinematic proxy           |
| `WE_identity_r`    | -     | Work-energy identity ratio; 1.0 is exact agreement between force work and CoM energy change |
| `foreaftAxis`      | -     | Which raw plate axis was resolved as anatomical fore-aft                          |

---

## steps_index.csv

One row per step (3842 rows, 14 columns). The compact step-level index; the full measures are in `perStep_long_multi.csv`.

| Column            | Units | Description                                                                          |
| :---------------- | :---- | :------------------------------------------------------------------------------------ |
| `boutID`          | -     | Trial identifier                                                                     |
| `stepIndex`       | -     | Step number within the trial                                                         |
| `startSample`     | -     | First sample of the step, indexing the per-trial time series                         |
| `endSample`       | -     | Last sample of the step                                                              |
| `startTime_s`     | s     | Time of step start                                                                   |
| `endTime_s`       | s     | Time of step end                                                                     |
| `gait`            | -     | `walk`, `groundedRun` or `aerialRun`; the classification the paper uses              |
| `gait4`           | -     | `walkGrounded`, `groundedRun`, `aerialRun` or `pendularRun`; as `gait`, with the aerial-pendular cell (the pendular run) split out of aerial running |
| `steadiness`      | -     | `steady`, `accelerating` or `decelerating`, per-stride grade inherited by both steps |
| `outlier`         | -     | Mis-cut cycle flag                                                                   |
| `qcPass`          | -     | Retention-gate flag                                                                  |
| `analysisSample`  | -     | 1 for the rows behind the paper                                                      |
| `exclusionReason` | -     | Why a row is not in the analysis sample                                              |
| `u`               | -     | Dimensionless speed, $v/\sqrt{gL_0}$                                                 |

---

## perStep_long_multi.csv

The tidy per-step analysis table: 3842 rows, 85 columns, one row per step. This and `perStride_long_multi.csv` are the tables the analysis and the summary statistics are computed from.

Columns ending `_n` are dimensionless (see Normalisation). Raw and dimensionless forms are both kept so a reader can work in either.

**Identity and grouping**

| Column | Description |
| :--- | :--- |
| `boutID`, `stepIndex` | Trial and step, the join key to `steps_index.csv` and the per-trial series |
| `bird`, `colourCode`, `study`, `dataset`, `dateCode`, `trialNo` | Session, cohort and collection labels |
| `mass_kg`, `L0_m`, `BW_N` | Per-session body mass, leg length and body weight, the normalisation scales |
| `gaitLabelHand` | Legacy hand label from the original collection; not used in the analysis |

**Classification and sample flags**

| Column | Description |
| :--- | :--- |
| `gait`, `gait4` | The classification the paper uses: aerial phase crossed with the direction of energy exchange |
| `gaitHodo`, `mechHodo` | The alternative classification by the sign of the CoM velocity-loop area, which the paper tests against `gait` |
| `steadiness`, `accClass`, `accClass_fa` | Steadiness class by the energy grade (analysis criterion), by each row's own energy grade, and by the fore-aft speed change |
| `analysisSample`, `exclusionReason`, `outlier` | Sample flags, as in `steps_index.csv` |

**Timing and duty**

| Column | Units | Description |
| :--- | :--- | :--- |
| `stepPeriod`, `stepPeriod_n` | s | Step duration |
| `stepFreq`, `stepFreq_n` | Hz | Step frequency |
| `contactTime`, `contactTime_n` | s | Marker-based stance duration |
| `flightTime`, `flightTime_n` | s | Aerial-phase duration |
| `hasFlight` | - | 1 where the summed vertical force reaches zero |
| `dutyFactor` | - | Per-limb duty factor from foot-marker contact; missing where marker contact could not be determined |
| `dutyFactor_gf` | - | Force-based estimate, grounded fraction divided by two, kept for cross-check |
| `groundedFrac` | - | Fraction of the step with any ground contact |
| `dutyLimb` | - | Which limb the duty factor refers to |

**Speed and distance**

| Column | Units | Description |
| :--- | :--- | :--- |
| `meanSpeed`, `meanSpeed_n` | m/s | Mean forward speed over the step. `meanSpeed_n` is the dimensionless speed $u$ |
| `speed_TD`, `speed_TD_n` | m/s | Forward speed at touchdown |
| `stepLength`, `stepLength_n` | m | Fore-aft CoM displacement over the step |
| `Froude` | - | Froude number, $u^2$ |
| `accelForeAft`, `accelForeAft_n`, `accelSign` | m/s² | Mean fore-aft acceleration and its sign |
| `vertExcursion`, `vertExcursion_n` | m | Vertical CoM excursion |

**Forces and impulse**

| Column | Units | Description |
| :--- | :--- | :--- |
| `peakVertForce_BW` | BW | Peak vertical ground reaction force |
| `peakForeAftForce_BW` | BW | Peak absolute fore-aft force |
| `peakResultantForce_BW` | BW | Peak resultant force magnitude |
| `impulseForeAft`, `impulseForeAft_n` | N s | Net fore-aft impulse over stance |

**Center-of-mass energy and work**

| Column | Units | Description |
| :--- | :--- | :--- |
| `dKE_step`, `dKE_step_n` | J | Change in forward kinetic energy $E_{kf} = \frac{1}{2}m(v_{ml}^2 + v_{fa}^2)$ |
| `dPE_step`, `dPE_step_n` | J | Change in vertical energy $E_v = mgz + \frac{1}{2}mv_{vert}^2$ |
| `dE_CoM`, `dE_CoM_n` | J | Net change in total CoM mechanical energy |
| `W_CoM_pos`, `W_CoM_neg`, `W_CoM_net` (and `_n`) | J | Positive, negative and net CoM work |
| `KE_range`, `PE_range` (and `_n`) | J | Within-step range of each energy component |
| `ampRatio_KEtoPE` | - | Ratio of the kinetic to the vertical energy amplitude |
| `WE_relresid`, `WE_r` | - | Work-energy identity relative residual and ratio, the reconstruction-quality measures |
| `fracEG`, `fracDE`, `fracDV` | - | Steadiness measures: energy grade $\Delta E_{CoM}/(mgL)$ (the analysis criterion), energy fraction $\Delta E_{CoM}/(mv^2)$, and fore-aft speed change $\lvert aT/v \rvert$ |
| `driftRMS_vert_mm` | mm | Trial-level vertical reconstruction drift, repeated on each step |

**Gait descriptors** (the six continuous descriptors of the gait space, plus the hodograph measure)

| Column | Units | Description |
| :--- | :--- | :--- |
| `recovery` | % | Pendular energy recovery (Cavagna, Heglund and Taylor 1977) |
| `congruity` | % | Percentage of the cycle over which $E_{kf}$ and $E_v$ change in the same direction (Ahn, Furrow and Biewener 2004) |
| `collisionAngle` | rad | Force- and velocity-weighted deviation of the ground reaction force from perpendicular to CoM velocity (Lee et al. 2011, 2013) |
| `collisionFraction` | - | Fraction of the cycle contributing to the collision angle |
| `CoTmech` | - | Dimensionless mechanical cost of transport |
| `hodoArea`, `hodoArea_n` | m²/s² | Signed area of the CoM velocity loop by the shoelace formula. Positive is counterclockwise (pendular), negative clockwise (bouncing). `hodoArea_n` is $A/(gL_0)$ |

**Virtual leg and leg stiffness**

| Column | Units | Description |
| :--- | :--- | :--- |
| `legLen_TD`, `legLen_TD_n` | m | CoM-to-foot distance at touchdown |
| `legAngle_TD` | deg | Virtual leg angle at touchdown, counterclockwise from horizontal with fore-aft positive in the direction of travel |
| `dLegLen`, `dLegAngle` | m, deg | Change in leg length and angle over stance; a negative `dLegAngle` is limb retraction |
| `legCompress_n` | - | Leg compression as a fraction of touchdown leg length |
| `kLeg_n` | BW/L0 | Dimensionless leg-spring stiffness $kL_0/(mg)$, from $k = F_{max}/\Delta L$ with compression recovered by inverting a sinusoidal vertical force profile over contact (Blum, Lipfert and Seyfarth 2009, method C) |

`legCompress_n` and `kLeg_n` are missing where the inversion returns a non-physical compression. Leg stiffness is defined only relative to the spring-mass template it comes from, and where two limbs share the net force in double support it is an effective two-limb stiffness rather than a single-limb property.

---

## perStride_long_multi.csv

The tidy per-stride table: 1678 rows, 77 columns, one row per stride (two steps). Column meanings match `perStep_long_multi.csv` except as noted.

| Column group | Columns | Notes |
| :--- | :--- | :--- |
| Identity | `boutID`, `strideIndex`, `bird`, `colourCode`, `study`, `dataset`, `dateCode`, `trialNo`, `mass_kg`, `L0_m`, `BW_N`, `gaitLabelHand`, `gaitCodeHand` | `strideIndex` joins to the per-step table via the two steps it spans |
| Sample flag | `analysisSample`, `outlier` | The stride table carries no `gait` column; join to the per-step table for gait |
| Timing | `stridePeriod`, `strideFreq`, `contactTime`, `flightTime`, `hasFlight`, `dutyFactor`, `dutyFactor_gf`, `groundedFrac`, `dutyLimb` (and `_n` forms) | Per stride rather than per step |
| Speed and distance | `meanSpeed`, `strideLength`, `Froude`, `accelForeAft`, `accelSign`, `vertExcursion` (and `_n` forms) | |
| Energy and work | `dKE`, `dPE`, `Wf_horiz`, `Wv_vert`, `Wtot_com`, `W_CoM_pos`, `W_CoM_neg`, `W_CoM_net`, `dE_CoM`, `KE_range`, `PE_range`, `ampRatio_KEtoPE`, `WE_relresid`, `WE_r`, `driftRMS_vert_mm` (and `_n` forms) | `Wf_horiz`, `Wv_vert` and `Wtot_com` are the summed positive increments of the horizontal, vertical and total energy that enter the recovery calculation |
| Steadiness | `fracEG`, `fracDE`, `fracDV`, `accClass`, `accClass_fa` | Steadiness is defined at stride grain; the per-step table inherits it |
| Descriptors | `recovery`, `congruity`, `collisionAngle`, `collisionFraction`, `CoTmech`, `hodoArea`, `hodoArea_n`, `gaitHodo`, `mechHodo` | |

Pendular recovery is $R = (W_f + W_v - W_{tot})/(W_f + W_v)$, each $W$ the sum of positive increments over the stride, and is defined per stride because the CoM returns to its state over a full stride.

---

## morphology_multi.csv

One row per study session (18 rows, 7 columns).

| Column     | Units | Description                                                              |
| :--------- | :---- | :------------------------------------------------------------------------ |
| `bird`     | -     | Study-session bird code                                                  |
| `study`    | -     | Study tag                                                                |
| `mass_kg`  | kg    | Body mass for that session                                               |
| `L0_m`     | m     | Leg length, the median CoM-to-foot distance at touchdown over the session's steps with no aerial phase |
| `Lnorm_m`  | m     | Geometric length scale, $0.20\,m^{1/3}$, reported for reference; it is not the normaliser used |
| `nSteps`   | -     | Steps detected for that session                                          |
| `nStrides` | -     | Strides detected for that session                                        |

---

## trialRoster_multi.csv

Source-file provenance, one row per trial (334 rows, 12 columns): `boutID`, `bird`, `colourCode`, `study`, `dateCode`, `trialNo`, `trialStem`, `bodyMass_kg`, `windowMode`, and the `forceFile`, `markerFile` and `logFile` paths of the raw recordings.

Paths are given relative to the root of the `LEVEL_collation` archive, for example
`OtherSpecies_byDate/Guinea fowl/Blum_Surface_2012-06/2012-06-25/gf_red_LE0cm_0006.txt`, so the roster is portable across machines. The pipeline locates that archive at run time, under the project folder or one level above it, and joins these paths to whichever it finds. `logFile` is empty for the trials that have none. The raw force and marker files are not part of this deposit; the roster records which source file each trial came from.

---

## mean_stepcycle_by_gait_speed

Ensemble mean step cycles for each gait by speed-bin group, in long format (10800 rows) and as a MATLAB struct. Six groups: walk, grounded run and aerial run, each split slow and fast.

Two things determine which cycles enter a group, which corresponds to the paper's mean-trace figure. 1) The cycles are the STEADY steps of the analysis sample (`steadiness` = `steady` and `analysisSample` = 1). 2) The slow/fast boundary is the gait's median Froude number over all of its classified steps, not the median of the steady subset, so the boundary does not move between steadiness classes; the two halves of a group can therefore be unequal in size.

| Column      | Units | Description                                                              |
| :---------- | :---- | :------------------------------------------------------------------------ |
| `gait`      | -     | `walk`, `groundedRun` or `aerialRun`                                     |
| `speedbin`  | -     | `slow` or `fast`, split at the gait's median Froude over all its steps    |
| `n`         | -     | Number of step cycles in the group mean                                  |
| `dur_mean_s`| s     | Mean cycle duration for the group                                        |
| `pct`       | %     | Percent of the step cycle, 0 to 100                                      |
| `time_s`    | s     | Real-time axis, `pct` scaled by `dur_mean_s`, preserving the mean duration |
| `channel`   | -     | Which of the 18 channels this row belongs to                             |
| `mean`      | varies| Group mean of that channel at that phase                                 |
| `sd`        | varies| Standard deviation across cycles                                         |

Channels: `F_vt_N`, `F_fa_N`, `F_ml_N`, `F_vt_BW`, `F_fa_BW` (forces, raw and body-weight normalised); `vfa_ms`, `vvert_ms` (CoM velocity); `com_fa_disp_m`, `com_vt_rel_m` (CoM displacement); `KE_J`, `PE_J`, `Etot_J`, `KE_n`, `PE_n` (energy, raw and normalised); `stanceFoot_fa_relCoM_m`, `stanceFoot_vt_relCoM_m`, `swingFoot_fa_relCoM_m`, `swingFoot_vt_relCoM_m` (foot position relative to the CoM).

---

## Coordinate Frame and Units

Every three-vector is ordered `[medio-lateral, fore-aft, vertical]`. The frame is bird-centered and right-handed: fore-aft is positive in the direction of travel for every trial, vertical is positive up, and medio-lateral completes the right-handed frame. Trials in which the bird crossed the runway in the opposite lab direction have their fore-aft and medio-lateral signs flipped together, so all trials are directly comparable without further reorientation.

Units are SI: force in newtons, position in meters, velocity in meters per second, time in seconds, mass in kilograms, energy in joules. Angles are in degrees except the collision angle, which is in radians because the equivalence with the dimensionless cost of transport holds in that form.

## Normalisation

Data are normalised into dimensionless quantities following dynamic similarity, using each session's body mass $m$, leg length $L_0$ and $g = 9.81$ m/s². Columns ending `_n` carry the dimensionless form.

| Quantity | Normalisation factor | Formula        |
| :------- | :------------------- | :------------- |
| Force    | Body weight          | $mg$           |
| Length   | Leg length           | $L_0$          |
| Velocity | Froude velocity      | $\sqrt{gL_0}$  |
| Time     | Pendular period      | $\sqrt{L_0/g}$ |
| Frequency| Inverse pendular period | $\sqrt{g/L_0}$ |
| Work, energy | Body weight times leg length | $mgL_0$ |
| Impulse  | Body weight times pendular period | $mg\sqrt{L_0/g}$ |
| Loop area| Gravity times leg length | $gL_0$      |

Dimensionless speed is $u = v/\sqrt{gL_0}$ and the Froude number is $u^2$.

---

## Processing Workflow

The analysis code is on GitHub at [https://github.com/NeuroMechLab/GuineaFowlGait-code](https://github.com/NeuroMechLab/GuineaFowlGait-code), tagged `v1.2.1`. `PIPELINE.md` there gives the full run order and maps every output to the manuscript element that uses it. The sequence that produces the tables in this package is:

### MATLAB phase

**Step 1: Build the trial roster.** `GaitSelMulti_BuildRoster.m` pairs the marker and force files of the three collections and seeds body mass, writing `trialRoster_multi.csv`.

**Step 2: Import, reconstruct, detect, measure.** `GaitSelMulti_BatchProcess.m` imports each trial, resolves the fore-aft axis, reconstructs the CoM, detects gait events at foot-marker touchdowns, and computes the per-step and per-stride measures.

**Step 3: Normalise and write the tidy tables.** `GaitSelMulti_ExportTidyCSV.m` derives per-session `L0`, adds the dimensionless and leg-stiffness columns, and flags outliers, writing `perStep_long.csv`, `perStride_long.csv` and `morphology.csv`.

**Step 4: Quality-control gate.** `GaitSelMulti_QCGate.m` applies the work-energy and drift criteria and the energy-grade steadiness criterion. This is the single definition of the retained set, computed once in MATLAB so that this package and the statistics agree.

**Step 5: Build this package.** `GaitSelMulti_ExportDryad.m` writes the per-trial time series, `trial_index.csv` and `steps_index.csv`. `GaitSelMulti_ExportCycleTraces.m` then writes the per-step cycle traces the figures read, and `GaitSel_ExportDryadMeanCycles.m` writes `mean_cycles`.

### R phase

**Step 6: Clean and fix the analysis sample.** `Rscript run_all.R` runs the numbered chain. `01_load.R` reads the tidy tables and builds the individual grouping; `02_clean.R` applies the MATLAB quality-control and steadiness flags; `03_step_descriptors.R` adds the collision angle and touchdown energy; `04_analysis_sample.R` fixes the analysis sample of 2580 steps and 944 strides and writes the label table.

**Step 7: Figures, tables and statistics.** The remaining numbered scripts produce them. `21_dryad_labels.R` writes this package's `perStep_long_multi.csv` and `perStride_long_multi.csv`, adding `gait`, `gait4`, `gaitHodo`, `steadiness`, `analysisSample` and `exclusionReason` to the MATLAB tidy tables.

Step 5 reads the label table that step 6 writes, so `steps_index.csv` carries the gait and sample columns after `GaitSelMulti_ExportDryad.m` is run once more following the R phase.

### Reconstruction method

The CoM is reconstructed by path-matched double integration of the ground reaction force over the whole on-plate bout: $F = ma$, with one set of integration constants per axis (initial position, initial velocity and an acceleration-baseline offset) fitted to minimise departure from the kinematic CoM proxy. The net force is low-pass filtered at 50 Hz with a zero-phase Butterworth filter. The analysis window for each trial is found automatically as the span where the bird is both tracked and on the plates. Because the two rig generations mount their plates differently relative to the motion-capture frame, the horizontal force axes are assigned to anatomical fore-aft and medio-lateral by matching the cumulative horizontal impulse to the kinematic travel direction, using a consensus across each collection date.

The proxy is the midpoint of the cranial and caudal back markers, which lies above and behind the whole-body CoM, so a fore-aft and a vertical offset are applied to it before any geometry is taken. The fore-aft offset is fitted by driving the net pitch impulse about the CoM to zero over stride-complete windows, within anatomical bounds set by the two back markers. That condition does not identify the vertical offset, because the fore-aft impulse of the summed force approaches zero over a whole stride, so the vertical offset is set from anatomical scaling as $-0.2625\,L_\mathrm{iso}$ with $L_\mathrm{iso} = 0.20\,m^{1/3}$. Offsets are fitted for each bird in each recording session, because marker placement was renewed at each session. `comProxy_*` in the per-trial series carries them, so leg length, leg angle and CoM height above the foot all refer to one CoM. Energy fluctuations are unaffected, because a constant shift cancels in a height change.

Reconstruction quality per study, as median work-energy identity $r$ and median vertical drift: RVC 0.98 and 7.8 mm, DropVsPothole 1.00 and 7.5 mm, Surface 0.96 and 4.9 mm.

---

## Code Dependencies

**MATLAB** R2019b or later, with the Statistics and Machine Learning, Signal Processing and Curve Fitting Toolboxes.

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

Daley, M. A. and Birn-Jeffery, A. (2026). Data and code from: A velocity-loop view of avian gait:
hodograph rotation sense characterizes the walking-to-running continuum in guinea fowl. Dryad,
doi:10.5061/dryad.7m0cfxqc4

Contact: For questions regarding this dataset or code, please contact the corresponding author Monica Daley ([madaley@uci.edu](mailto:madaley@uci.edu)).
