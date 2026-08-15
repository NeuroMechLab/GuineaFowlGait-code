%% GaitSelMulti_RunFullWorkflow.m  —  MASTER SCRIPT, GF GAIT-SELECTION PIPELINE
% =========================================================================
% Guinea fowl gait-selection analysis over the three level-running collections (RVC 2008-2011 +
% Blum Surface 2012-06 + Blum DropVsPothole 2012-02). START HERE to regenerate EVERY figure,
% table, and the Data Dryad package from the raw force/marker files.
%
% The pipeline runs as ONE MATLAB phase (STEPS 1-6) followed by a SINGLE hand-off to R
% (STEP 7). Each %% block is independently re-runnable. Requires MATLAB R2019b+ (Statistics,
% Signal Processing, Curve Fitting) and R with the packages run_all.R loads. The raw
% force/marker files in the LEVEL_collation archive are needed only if the import cache is
% rebuilt in STEP 2. projectPaths locates that archive under the project root or one level
% above it; it currently sits one level above, shared with the parallel ostrich projects.
%
% DEPENDENCY CHAIN (why this order): the QC gate (qcPass) is computed in MATLAB (STEP 4,
% GaitSelMulti_QCGate), so the Dryad export can embed it (STEP 5) without waiting on R. The
% cycle-trace, mean-cycle and mean-trace exporters read the Dryad per-trial files, so they
% follow the Dryad export. Once all MATLAB products are final, R runs once (STEP 7): 02_clean.R
% READS the MATLAB qcPass/steadiness (single source of truth) and every figure/table is built
% from final inputs. 
% =========================================================================
P = projectPaths(); addpath(P.code); addpath(P.helpers);
rdir = fullfile(P.root, '_RAnalysis');
fprintf('Project root: %s\n', P.root);
% Resolve Rscript for the single R hand-off. MATLAB's system() uses a non-login shell whose
% PATH may lack Homebrew, so fall back to common install locations; set RSCRIPT by hand if yours
% is elsewhere.
RSCRIPT = 'Rscript';
if system('command -v Rscript >/dev/null 2>&1') ~= 0
    cands = {'/opt/homebrew/bin/Rscript','/usr/local/bin/Rscript','/usr/bin/Rscript'};
    hit = cands(cellfun(@isfile, cands));
    if ~isempty(hit), RSCRIPT = hit{1}; end
end

%% STEP 1 (MATLAB) — Build the bout roster (REQUIRED when the data set changes)
% Walks the level-trial folders of the three collections, study-tags the birds, seeds mass
% from AllBirds_metadata.csv. Output: _RAnalysis/data/trialRoster_multi.csv.
GaitSelMulti_BuildRoster();

%% STEP 2 (MATLAB) — Import, reconstruct, detect, measure every trial (REQUIRED)
% Import (cached) -> per-(study,date) consensus fore-aft axis -> whole-bout path-matched CoM
% -> gait events (foot-marker touchdowns) -> per-step/stride measures + objective gait class.
% Output: SavedResults/perStepMeasures_multi.mat + importCache_multi.mat, perStep_raw_multi.csv,
% perStride_raw_multi.csv.
%   Rebuild the import cache from the raw files with:  GaitSelMulti_BatchProcess(struct('useCache',false));
GaitSelMulti_BatchProcess();

%% STEP 3 (MATLAB) — Normalise + write the tidy hand-off tables (REQUIRED)
% Per-bird L0, dimensionless columns, and the speed-conditional outlier flag. Writes the
% *_multi.csv tables AND promotes them to the canonical names the R pipeline reads
% (perStep_long.csv, perStride_long.csv, morphology.csv). Also adds the per-step leg-stiffness
% columns (legCompress_n, kLeg_n) via addLegStiffness, which needs the per-bird L0 and so runs
% here rather than in the per-step measures.
GaitSelMulti_ExportTidyCSV();

%% STEP 4 (MATLAB) — QC gate + steadiness (REQUIRED; single source of truth)
% Work-energy residual + duration-normalized drift -> step_qcpass.csv / stride_qcpass.csv;
% net-CoM-energy steadiness -> step_steadiness.csv. The Dryad export and 02_clean.R both read
% these, so the analyzed set has ONE definition.
GaitSelMulti_QCGate();

%% STEP 5 (MATLAB) — Data Dryad package + cycle-trace and mean-cycle products (REQUIRED)
% (a) Per-trial time series + trial_index + steps_index (with qcPass) + cleaned foot tracks;
% (b) per-step cycle traces (cycleTracesStep.csv), read by the R step-descriptor stage and by
%     the hodograph, energy-exchange and mean-trace figures; (c) gait x speed mean step cycles
%     for the Dryad package. Both (b) and (c) read the per-trial series written by (a).
% (a) also reads _RAnalysis/data/step_labels_analysis.csv, which STEP 7 writes, so steps_index
% carries the gait and analysis-sample columns once this step is re-run after STEP 7.
%
% RE-RUNNING THIS STEP AFTER STEP 7 REVERTS TWO PACKAGE FILES. (a) copies the tidy tables into
% the package as perStep_long_multi.csv and perStride_long_multi.csv, which drops the label
% columns that R's 21_dryad_labels.R adds (gait, gait4, steadiness, analysisSample,
% exclusionReason, and the gaitObjective -> gaitHodo rename). So the order is STEP 5, STEP 7,
% STEP 5 again for steps_index, then 21_dryad_labels.R once more. That script is idempotent: it
% builds from data/, not by mutating the package copies. Check gait4 in
% GaitSel_DryadPackage_AllGF/perStep_long_multi.csv afterwards; if the column is absent, the
% relabel has not been re-run.
%
% (c) reads _RAnalysis/output/speedbin_thresholds_R.csv, also written by STEP 7, so the mean
% cycles are binned on the same slow/fast threshold as the paper's figures and will error rather
% than guess if the R pipeline has not run.
GaitSelMulti_ExportDryad();
GaitSelMulti_ExportCycleTraces();
GaitSel_ExportDryadMeanCycles(struct('root', fullfile(P.root, 'GaitSel_DryadPackage_AllGF')));

%% STEP 6 (MATLAB) — Multi-rig QC figures + per-trial manifest (validation, OPTIONAL)
% Raw-vs-filtered force, reconstruction-vs-kinematics, event overlays, step-duration
% distributions per rig, and manifest.csv. Written to _Internal/pipeline_qc/multi/. These are
% human-review diagnostics, not a paper product, so the script lives in _Internal and the call
% is skipped when it is absent (as in the publication bundle).
qcDir = fullfile(P.root, '_Internal', 'pipeline_qc');
if isfolder(qcDir)
    addpath(qcDir);
    GaitSel_QCPlots_Multi();
else
    fprintf('STEP 6 skipped: _Internal/pipeline_qc absent (internal QC figures).\n');
end

%% STEP 7 (R) — Single hand-off: all figures, tables and statistics (REQUIRED)
% One pass of run_all.R now that every input is final (tidy tables, qcPass/steadiness,
% cycleTracesStep.csv). The status is CHECKED, so a missing Rscript (for
% example when launched from the MATLAB app with a minimal PATH) errors here rather than
% silently skipping. If it errors, run it by hand from _RAnalysis:  Rscript run_all.R
status = system(sprintf('cd "%s" && "%s" run_all.R', rdir, RSCRIPT));
assert(status == 0, 'GaitSelMulti:Rfailed', ...
    ['R hand-off failed (status %d). Ensure Rscript is callable (tried "%s"); ' ...
     'then run by hand from _RAnalysis:  Rscript run_all.R'], status, RSCRIPT);

fprintf('\nWorkflow complete. Figures/tables: _RAnalysis/output/ ; Dryad package: GaitSel_DryadPackage_AllGF/\n');
