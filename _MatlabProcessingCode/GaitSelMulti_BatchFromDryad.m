function GaitSelMulti_BatchFromDryad(opts)
%GAITSELMULTI_BATCHFROMDRYAD  Re-run the measurement chain from the deposited per-trial series.
%
%   GaitSelMulti_BatchFromDryad()
%   GaitSelMulti_BatchFromDryad(struct('root', <package dir>, 'limit', 20))
%
%   The deposited per-trial files carry everything the measurement kernels read: the net ground
%   reaction force, the time base, the kinematic centre-of-mass proxy the reconstruction is
%   path-matched to, both foot marker tracks, and body mass, leg length and g. So the deposit is
%   the analysis input, and the reconstruction, the gait-event detection and the per-step
%   reduction all run again from it without the raw recordings.
%
%   Each trial is read, rebuilt into the bout struct the kernels expect, and passed through
%   GaitSel_ReconstructCoM, GaitSel_DetectGaitEvents and GaitSel_PerStepStrideMeasures. The
%   result is written as SavedResults/perStepMeasures_multi.mat, which is what
%   GaitSelMulti_ExportTidyCSV reads, so every step after this one is unchanged.
%
%   NOTHING IS WRITTEN INTO THE DEPOSIT. The downloaded package is opened read-only and the
%   outputs go to rerun/ and SavedResults/, so a second run starts from the archived state
%   rather than from the first run's results.
%
%   rerun/reconstruction_check.csv compares this run against the values the deposit carries for
%   the same trial: the vertical centre-of-mass path, the drift, the work-energy identity and the
%   number of steps cut. That is the check that the re-run reproduces the published analysis.
%
%   ONE COLUMN CANNOT BE REPRODUCED. dutyLimb is estimated from the per-plate vertical force,
%   which the deposit does not carry, so it comes back empty. The paper's duty factor is the
%   marker-based one and is unaffected.
%
%   See also: GaitSelMulti_ExportTidyCSV, GaitSelMulti_QCGate, GaitSel_ReconstructCoM.

    if nargin < 1, opts = struct(); end
    P     = projectPaths();
    root  = getOpt(opts, 'root', P.dryad);
    tsDir = fullfile(root, 'per_trial_timeseries');
    limit = getOpt(opts, 'limit', Inf);

    if ~isfolder(tsDir)
        error('GaitSelMulti_BatchFromDryad:noDeposit', ...
              ['per_trial_timeseries not found at\n    %s\n' ...
               'Download the Dryad package (doi:10.5061/dryad.7m0cfxqc4), unzip ' ...
               'per_trial_timeseries.zip inside it, and put the package at\n    %s\n' ...
               'or pass its location as opts.root.'], tsDir, root);
    end

    files = dir(fullfile(tsDir, '*.mat'));
    if isempty(files)
        error('GaitSelMulti_BatchFromDryad:emptyDeposit', ...
              'no .mat trials in %s; unzip per_trial_timeseries.zip first', tsDir);
    end
    n = min(numel(files), limit);
    fprintf('re-running %d of %d deposited trials from\n  %s\n', n, numel(files), tsDir);

    massBird = loadSessionMass(root);

    stepRows = {}; strideRows = {}; chk = {};
    faAxisUsed = strings(n, 1);
    for i = 1:n
        f = fullfile(tsDir, files(i).name);
        try
            B  = boutFromDeposit(f);
            dep = depositedRefs(f);
            B  = GaitSel_ReconstructCoM(B);
            E  = GaitSel_DetectGaitEvents(B);
            [sr, dr] = GaitSel_PerStepStrideMeasures(B, E);
            sr = tagRows(sr, B); dr = tagRows(dr, B);
            stepRows   = [stepRows, sr];    %#ok<AGROW>
            strideRows = [strideRows, dr];  %#ok<AGROW>
            faAxisUsed(i) = "deposited";    % the frame is already anatomical in the deposit
            chk{end+1} = compareToDeposit(B, E, dep); %#ok<AGROW>
        catch ME
            warning('%s: %s', files(i).name, ME.message);
        end
        if mod(i, 40) == 0, fprintf('  measured %d/%d\n', i, n); end
    end

    T = struct2table([stepRows{:}]);
    S = struct2table([strideRows{:}]);
    if ~isfolder(P.saved), mkdir(P.saved); end
    if ~isfolder(P.rerun), mkdir(P.rerun); end
    save(fullfile(P.saved, 'perStepMeasures_multi.mat'), 'S', 'T', 'massBird', 'faAxisUsed');
    writetable(T, fullfile(P.rData, 'perStep_raw_multi.csv'));
    writetable(S, fullfile(P.rData, 'perStride_raw_multi.csv'));

    C = struct2table([chk{:}]);
    writetable(C, fullfile(P.rerun, 'reconstruction_check.csv'));
    fprintf(['\n%d steps, %d strides over %d trials\n' ...
             'agreement with the deposit, worst trial:\n' ...
             '  vertical CoM path   %.3f mm\n' ...
             '  drift RMS           %.3f mm\n' ...
             '  work-energy ratio   %.5f\n' ...
             '  trials whose step count differs  %d of %d\n' ...
             'per-trial detail: %s\n'], ...
            height(T), height(S), height(C), ...
            max(C.maxAbsComVert_mm), max(abs(C.driftRMS_diff_mm)), ...
            max(abs(C.WE_r_diff)), sum(C.nSteps_rerun ~= C.nSteps_deposit), height(C), ...
            fullfile(P.rerun, 'reconstruction_check.csv'));
end

% ------------------------------------------------------------------------ helpers
function B = boutFromDeposit(f)
%BOUTFROMDEPOSIT  The bout struct the kernels read, from one deposited trial file.
    d  = load(f); tr = d.trial;
    B  = struct();
    B.boutID   = char(tr.boutID);
    B.bird     = string(tr.bird);
    B.study    = string(tr.study);
    B.dateCode = tr.dateCode;
    B.colourCode = string(getField(tr, 'colourCode', ""));
    B.mass     = tr.mass_kg;
    B.g        = tr.g;
    B.fHz      = tr.fs_force_Hz;
    B.kinHz    = getField(tr, 'fs_kin_Hz', NaN);
    B.nPlates  = getField(tr, 'nPlates', NaN);
    B.time     = tr.time_s(:);
    B.force    = asNby3(tr.GRF_N);
    B.comProxy = asNby3(tr.comProxy_m);
    B.footR    = asNby3(tr.footR_m);
    B.footL    = asNby3(tr.footL_m);
    % legacy hand labels are not part of the deposit and are not used by the analysis
    B.gaitLabel = ""; B.gaitCode = NaN; B.trialNo = NaN;
end

function d = depositedRefs(f)
%DEPOSITEDREFS  The published results for this trial, to compare the re-run against.
    s = load(f); tr = s.trial;
    d.com      = asNby3(tr.com_m);
    d.driftRMS = tr.driftRMS_vert_mm;
    d.WE_r     = tr.WE_identity_r;
    d.nSteps   = numel(unique(tr.stepIndexPerSample(isfinite(tr.stepIndexPerSample) & ...
                                                   tr.stepIndexPerSample > 0)));
end

function c = compareToDeposit(B, E, dep)
%COMPARETODEPOSIT  One row of the agreement report.
    c.boutID = string(B.boutID);
    nv = min(size(B.com,1), size(dep.com,1));
    c.maxAbsComVert_mm = max(abs(B.com(1:nv,3) - dep.com(1:nv,3))) * 1000;
    % the kernel returns metres; the deposit stores millimetres
    rerunDrift_mm = B.driftRMS(3) * 1000;
    c.driftRMS_rerun_mm    = rerunDrift_mm;
    c.driftRMS_deposit_mm  = dep.driftRMS;
    c.driftRMS_diff_mm     = rerunDrift_mm - dep.driftRMS;
    c.WE_r_rerun    = B.WEcheck.r;
    c.WE_r_deposit  = dep.WE_r;
    c.WE_r_diff     = B.WEcheck.r - dep.WE_r;
    c.nSteps_rerun   = max(numel(E.stepIdx) - 1, 0);
    c.nSteps_deposit = dep.nSteps;
end

function M = loadSessionMass(root)
%LOADSESSIONMASS  bird|YYYY-MM -> body mass, from the deposit's morphology table.
    p = fullfile(root, 'morphology_multi.csv');
    if ~isfile(p)
        error('GaitSelMulti_BatchFromDryad:noMorphology', ...
              'morphology_multi.csv not found in %s', root);
    end
    Tm = readtable(p, 'TextType', 'string');
    M  = containers.Map(cellstr(Tm.session), num2cell(Tm.mass_kg));
end

function rows = tagRows(rows, B)
    for j = 1:numel(rows)
        rows{j}.study      = B.study;
        rows{j}.dataset    = "archive_multi";
        rows{j}.colourCode = B.colourCode;
        rows{j}.dateCode   = B.dateCode;
    end
end

function X = asNby3(X)
%ASNBY3  Three-column orientation, whichever way the array was stored.
    if size(X,2) ~= 3 && size(X,1) == 3, X = X.'; end
end

function v = getField(s, n, d)
    if isstruct(s) && isfield(s, n) && ~isempty(s.(n)), v = s.(n); else, v = d; end
end

function v = getOpt(o, n, d)
    if isstruct(o) && isfield(o, n) && ~isempty(o.(n)), v = o.(n); else, v = d; end
end
