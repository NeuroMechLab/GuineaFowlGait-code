function [S, T] = GaitSelMulti_BatchProcess(opts)
%GAITSELMULTI_BATCHPROCESS  Process every level trial to per-step / per-stride
%   measures. Reads trialRoster_multi.csv, imports each bout with
%   GaitSelMulti_ImportBout, applies the per-(study,date) consensus fore-aft force
%   axis, reconstructs the CoM, detects gait events, reduces to per-step and
%   per-stride measures, tags each row with study/dataset, and writes:
%       SavedResults/perStepMeasures_multi.mat   (S, T, massBird, faAxisUsed)
%       _RAnalysis/data/perStep_raw_multi.csv
%       _RAnalysis/data/perStride_raw_multi.csv
%
%   Reports per-study reconstruction QC (work-energy identity r, vertical drift)
%   and the auto-calibrated fore-aft force axis, so within-rig consistency can be
%   checked.
%
%   See also: GaitSelMulti_BuildRoster, GaitSelMulti_ImportBout, GaitSelMulti_ExportTidyCSV.

    if nargin < 1, opts = struct(); end
    useCache = getOpt(opts,'useCache',true);
    g = 9.81;
    P = projectPaths(); addpath(P.helpers); addpath(P.code);
    % Delimiter is pinned: the roster's path fields contain spaces, and readtable's
    % automatic delimiter detection can prefer the space over the comma.
    R = readtable(fullfile(P.rData,'trialRoster_multi.csv'), 'TextType','string', ...
                  'Delimiter',',');
    nB = height(R);

    % ---- import (cached) + per-bout mean vertical force + force cal -----
    cacheFile = fullfile(P.saved,'importCache_multi.mat');
    if useCache && exist(cacheFile,'file')
        load(cacheFile,'bouts','meanFz','calCx','calCy'); fprintf('Loaded multi import cache (%d bouts).\n',numel(bouts));
    else
        % ---- per-bird-session centre-of-mass offsets, resolved before import ----
        % Fitted by AvianGaitCore/matlab/runFitComOffsets.m from the net-zero pitch-impulse
        % condition over whole strides, under anatomical bounds, with the vertical from the
        % anatomical scaling because the condition cannot identify it. One row per BIRD-SESSION,
        % keyed birdCode|YYYY-MM: the markers are re-applied at each session, and the two 2012
        % sessions carry back markers about 4 cm apart. Documented in
        % ../_Documentation/PITCH_CORRECTION_METHOD.md.
        OFF = loadComOffsets(P);
        R.comOffsetFA   = zeros(nB,1);
        R.comOffsetVert = zeros(nB,1);
        R.comOffsetSrc  = strings(nB,1);
        for i = 1:nB
            k = offsetKey(R.bird(i), R.dateCode(i));
            if isKey(OFF,k)
                o = OFF(k);
                R.comOffsetFA(i) = o.fa; R.comOffsetVert(i) = o.vert; R.comOffsetSrc(i) = o.src;
            else
                R.comOffsetSrc(i) = "none";
                warning('GaitSelMulti:noOffset','no CoM offset for %s', k);
            end
        end
        fprintf('CoM offsets applied to %d/%d trials\n', sum(R.comOffsetVert~=0), nB);

        bouts = cell(nB,1); meanFz = nan(nB,1); calCx = nan(nB,1); calCy = nan(nB,1);
        for i = 1:nB
            try
                row = R(i,:);
                row.comOffset = [0, R.comOffsetFA(i), R.comOffsetVert(i)];   % [ml fa vert]
                B = GaitSelMulti_ImportBout(row);
                bouts{i} = B; meanFz(i) = mean(B.force(:,3),'omitnan');
                calCx(i) = B.forceCal.cx; calCy(i) = B.forceCal.cy;
            catch ME
                warning('import failed for %s: %s', R.boutID(i), ME.message);
            end
            if mod(i,25)==0, fprintf('  imported %d/%d\n', i, nB); end
        end
        save(cacheFile,'bouts','meanFz','calCx','calCy','-v7.3');
    end

    % ---- per-(study,date) consensus fore-aft force axis -----------------
    % One collection date = one rig setup, so the horizontal plate axis carrying
    % fore-aft is fixed within a (study,date); only the travel DIRECTION (sign)
    % varies per trial. Per-trial correlations are individually weak, so pick the
    % axis by aggregate evidence (larger summed |correlation|) across the group,
    % then take the sign per trial. Robust to the occasional ambiguous trial.
    grpKey = R.study + "|" + string(R.dateCode);
    [ug,~,gi] = unique(grpKey);
    faUseX = false(numel(ug),1);
    for k = 1:numel(ug)
        sel = gi==k;
        faUseX(k) = sum(abs(calCx(sel)),'omitnan') >= sum(abs(calCy(sel)),'omitnan');
    end

    % ---- mass from force, per BIRD-SESSION ------------------------------
    % Keyed birdCode|YYYY-MM, the same session key the CoM offsets are fitted on, so the mass
    % that scales a trial's reconstruction and the offset applied to it come from the same
    % recording session.
    sessKey  = R.bird + "|" + sessionMonth(R.dateCode);
    massBird = containers.Map('KeyType','char','ValueType','double');
    us = unique(sessKey);
    for i = 1:numel(us)
        sel = sessKey==us(i) & isfinite(meanFz);
        if any(sel); massBird(char(us(i))) = median(meanFz(sel))/g; end
    end

    % ---- reconstruct, detect, reduce to measures ------------------------
    stepRows = {}; strideRows = {}; faAxisUsed = strings(nB,1);
    for i = 1:nB
        B = bouts{i}; if isempty(B), continue; end
        sk = char(sessKey(i));
        if ~isKey(massBird,sk); continue; end
        B.mass = massBird(sk);
        % override the per-trial fore-aft axis with the (study,date) consensus
        Fx = B.forceRawHoriz(:,1); Fy = B.forceRawHoriz(:,2); Fz = B.forceVert;
        if faUseX(gi(i))
            B.force = [Fy, sgn(calCx(i))*Fx, Fz]; faAxisUsed(i) = "X";
        else
            B.force = [Fx, sgn(calCy(i))*Fy, Fz]; faAxisUsed(i) = "Y";
        end
        try
            B = GaitSel_ReconstructCoM(B);
            E = GaitSel_DetectGaitEvents(B);
            [sr, dr] = GaitSel_PerStepStrideMeasures(B, E);
            sr = tagRows(sr, R(i,:)); dr = tagRows(dr, R(i,:));
            stepRows   = [stepRows, sr];       %#ok<AGROW>
            strideRows = [strideRows, dr];     %#ok<AGROW>
        catch ME
            warning('reconstruct/measure failed for %s: %s', B.boutID, ME.message);
        end
        if mod(i,40)==0, fprintf('  measured %d/%d\n', i, nB); end
    end

    T = struct2table([stepRows{:}]);
    S = struct2table([strideRows{:}]);
    if ~exist(P.saved,'dir'); mkdir(P.saved); end
    save(fullfile(P.saved,'perStepMeasures_multi.mat'),'S','T','massBird','faAxisUsed');
    writetable(T, fullfile(P.rData,'perStep_raw_multi.csv'));
    writetable(S, fullfile(P.rData,'perStride_raw_multi.csv'));

    fprintf('\nBatchProcessMulti: %d steps, %d strides across %d bouts.\n', height(T), height(S), nB);
    % per-study QC report
    fprintf('\nPer-study QC (reconstruction + force-axis calibration):\n');
    studies = unique(R.study);
    for k = 1:numel(studies)
        bi = R.study==studies(k);
        ax = faAxisUsed(bi); ax = ax(ax~="");
        we = S.WE_r(ismember(S.boutID, R.boutID(bi)));
        dr = S.driftRMS_vert_mm(ismember(S.boutID, R.boutID(bi)));
        [ua,~,uc] = unique(ax); cnt = accumarray(uc,1);
        axstr = strjoin(arrayfun(@(a,c) sprintf('%s:%d',a,c), ua, cnt, 'uni',0), ' ');
        fprintf('  %-28s WE r med %.3f ; drift med %.1f mm ; fore-aft axis {%s}\n', ...
            studies(k), median(we,'omitnan'), median(dr,'omitnan'), axstr);
    end
end

function s = sgn(x)
    if ~isfinite(x) || x >= 0, s = 1; else, s = -1; end
end

% -------------------------------------------------------------------- helpers
function rows = tagRows(rows, rr)
    st = rr.study; ds = "archive_multi"; cc = rr.colourCode; sd = rr.dateCode;
    for j = 1:numel(rows)
        rows{j}.study = st; rows{j}.dataset = ds; rows{j}.colourCode = cc; rows{j}.dateCode = sd;
    end
end

function v = getOpt(o,n,d)
    if isstruct(o)&&isfield(o,n)&&~isempty(o.(n)); v=o.(n); else; v=d; end
end

function OFF = loadComOffsets(P)
%LOADCOMOFFSETS  Per-bird-session fore-aft and vertical CoM offsets from the archive.
%   Keyed birdCode|YYYY-MM, the fitted table's own fitKey, so a bird recorded in more than one
%   session takes the offset belonging to that session's marker placement.
    OFF = containers.Map('KeyType','char','ValueType','any');
    f = fullfile(P.collation,'_metadata','CoM_offsets_fitted.csv');
    if ~isfile(f)
        warning('GaitSelMulti:noOffsetFile','no %s; running with no CoM correction', f); return
    end
    opts = detectImportOptions(f,'Delimiter',',');
    opts = setvartype(opts,'string');
    opts.VariableNamingRule = 'preserve';
    T = readtable(f, opts);
    T = T(T.species == "Guinea fowl", :);
    for i = 1:height(T)
        OFF(char(strtrim(T.fitKey(i)))) = struct( ...
            'fa',   str2double(T.fa_offset_m(i)), ...
            'vert', str2double(T.vert_offset_m(i)), ...
            'src',  strtrim(T.source(i)));
    end
end

function k = offsetKey(bird, dateCode)
%OFFSETKEY  The fitted table's fitKey for one trial: bird code and session month.
    k = char(string(bird) + "|" + sessionMonth(dateCode));
end

function m = sessionMonth(dateCode)
%SESSIONMONTH  YYYY-MM from a YYYYMMDD date code, elementwise.
%   The session is written to month resolution because that is the scope the archive's offset
%   tables are keyed on. Month separates every guinea fowl session: no bird was recorded on two
%   dates within one month.
    d = string(dateCode);
    t = regexprep(d, '\D', '');
    assert(all(strlength(t) >= 6), 'GaitSelMulti:badDateCode', ...
           'a dateCode has no YYYYMM to key the session on');
    m = extractBetween(t,1,4) + "-" + extractBetween(t,5,6);
end
