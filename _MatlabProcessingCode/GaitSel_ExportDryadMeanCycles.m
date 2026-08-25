function GaitSel_ExportDryadMeanCycles(opts)
%GAITSEL_EXPORTDRYADMEANCYCLES  Data Dryad package, Phase 2: average step-cycle
%   trajectories by gait x speed bin. For each bin (walk / grounded run / aerial
%   run, each split slow/fast at the gait's median dimensionless speed), the mean
%   and s.d. over the 0-100% step cycle of: ground reaction force (3 axes, N and
%   body weights), CoM velocity (fore-aft, vertical), vertical CoM position (about
%   the step mean), fore-aft CoM displacement, each foot relative to the CoM
%   (fore-aft, vertical), and CoM energies (kinetic, potential, total; in J and
%   dimensionless /(m g L0)). Written as mean_cycles.mat and mean_cycles.csv.
%
%   Built from the per-trial files of the data package. Uses the
%   STEADY steps only (the clean cyclic set); each step is resampled to 100 points.
%   Axis order of any 3-vector is [medio-lateral, fore-aft, vertical].
    if nargin < 1, opts = struct(); end
    P = projectPaths();
    root  = getOpt(opts, 'root', fullfile(P.root, 'GaitSel_DryadPackage_AllGF'));
    tsDir = fullfile(root, 'per_trial_timeseries');
    outDir = fullfile(root, 'mean_stepcycle_by_gait_speed');
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    g = 9.81; NP = 100; gaits = {'walk', 'groundedRun', 'aerialRun'};

    files = dir(fullfile(tsDir, '*.mat'));
    % Foot channels are the STANCE (contact) foot and the SWING foot, not right/left:
    % step cycles pool successive (alternating-foot) steps under left-right symmetry, so
    % a fixed right/left average would mix stance and swing. Per step the stance foot is
    % the one in ground contact for the majority of the step, by the foot-marker contact
    % criterion used in the pipeline (foot height below its 5th percentile + 3 cm).
    chans = {'F_ml_N','F_fa_N','F_vt_N','F_fa_BW','F_vt_BW','vfa_ms','vvert_ms', ...
             'com_vt_rel_m','com_fa_disp_m','stanceFoot_fa_relCoM_m','stanceFoot_vt_relCoM_m', ...
             'swingFoot_fa_relCoM_m','swingFoot_vt_relCoM_m','KE_J','PE_J','Etot_J','KE_n','PE_n'};
    % uBin is the step's dimensionless speed AS THE TABLES CARRY IT, kept separate from the `u`
    % recomputed below. The recomputed one differs in the last digits, which moved one aerial-run
    % step across the slow/fast threshold and made this export disagree with the figures by one.
    S = struct('gait', {}, 'u', {}, 'uBin', {}, 'dur', {}, 'cyc', {});   % per-step store
    for f = 1:numel(files)
        L = load(fullfile(tsDir, files(f).name), 'trial'); B = L.trial;
        m = B.mass_kg; L0 = B.L0_m; BW = m * g; if ~(L0 > 0), continue; end
        st = B.steps; t = B.time_s; F = B.GRF_N; com = B.com_m; vel = B.comVel_ms;
        fR = B.footR_m; fL = B.footL_m; n = numel(t);
        % per-foot ground contact (foot height below its 5th percentile + 3 cm)
        cR = fR(:,3) <= prctile(fR(:,3), 5) + 0.03;
        cL = fL(:,3) <= prctile(fL(:,3), 5) + 0.03;
        % travel-direction sign for the fore-aft axis: bouts cross the runway both
        % ways, so pooling steps requires flipping fore-aft to the direction of travel
        % (otherwise opposite-direction bouts cancel). Matches the trace exporters.
        td = sign(com(end,2) - com(1,2)); if td == 0, td = 1; end
        for k = 1:height(st)
            if ~ismember(char(st.gait(k)), gaits), continue; end
            if ~strcmp(char(st.steadiness(k)), 'steady'), continue; end
            % The ANALYSIS SAMPLE, not merely the QC gate. The paper's sample is the QC gate plus
            % every gait-space descriptor defined (04_analysis_sample.R), so filtering on it makes
            % the package's mean cycles describe the set the paper's tables report.
            if ismember('analysisSample', st.Properties.VariableNames)
                if st.analysisSample(k) ~= 1, continue; end
            elseif ismember('qcPass', st.Properties.VariableNames) && st.qcPass(k) ~= 1
                continue
            end
            a = st.startSample(k); b = min(st.endSample(k), n); if b - a < 5, continue; end
            idx = a:b; np = numel(idx); ph = linspace(0, 100, NP);
            rs = @(xseg) interp1(linspace(0, 100, np), xseg(:)', ph, 'linear');   % xseg is the step segment
            comfa = com(idx, 2); comvt = com(idx, 3);
            dur = t(b) - t(a); u = abs(comfa(end) - comfa(1)) / max(dur, eps) / sqrt(g * L0);
            c = struct();
            c.F_ml_N = rs(F(idx,1)); c.F_fa_N = td * rs(F(idx,2)); c.F_vt_N = rs(F(idx,3));
            c.F_fa_BW = c.F_fa_N / BW; c.F_vt_BW = c.F_vt_N / BW;
            c.vfa_ms = td * rs(vel(idx,2)); c.vvert_ms = rs(vel(idx,3));
            c.com_vt_rel_m = rs(comvt - mean(comvt));
            c.com_fa_disp_m = td * rs(comfa - comfa(1));
            if sum(cR(idx)) >= sum(cL(idx)), stF = fR; swF = fL; else, stF = fL; swF = fR; end
            c.stanceFoot_fa_relCoM_m = td * rs(stF(idx,2) - com(idx,2)); c.stanceFoot_vt_relCoM_m = rs(stF(idx,3) - com(idx,3));
            c.swingFoot_fa_relCoM_m  = td * rs(swF(idx,2) - com(idx,2)); c.swingFoot_vt_relCoM_m  = rs(swF(idx,3) - com(idx,3));
            KEseg = 0.5 * m * (vel(idx,1).^2 + vel(idx,2).^2 + vel(idx,3).^2);
            PEseg = m * g * (comvt - mean(comvt));
            c.KE_J = rs(KEseg - mean(KEseg)); c.PE_J = rs(PEseg); c.Etot_J = c.KE_J + c.PE_J;
            c.KE_n = c.KE_J / (m * g * L0); c.PE_n = c.PE_J / (m * g * L0);
            uBin = u;   % fall back to the recomputed speed only if the table lacks the column
            if ismember('u', st.Properties.VariableNames) && isfinite(st.u(k)), uBin = st.u(k); end
            S(end+1) = struct('gait', char(st.gait(k)), 'u', u, 'uBin', uBin, ...
                              'dur', dur, 'cyc', c); %#ok<AGROW>
        end
    end

    % Slow/fast threshold: the gait's median Froude over ALL classified steps, read from
    % _RAnalysis/output/speedbin_thresholds_R.csv, which 05_mean_traces.R writes and the paper's
    % figures use. Fixing the threshold over the full set holds it steady across the steadiness
    % classes, and reading it from that one file makes the package's bins the figure's bins.
    % Compared in FROUDE space (u^2 against medianFroude), the comparison R makes, rather than
    % thresholding u at sqrt(medianFroude): the square root rounds, and a step sitting on the
    % threshold then falls the other way.
    med = containers.Map();
    thrFile = fullfile(P.root, '_RAnalysis', 'output', 'speedbin_thresholds_R.csv');
    if isfile(thrFile)
        T = readtable(thrFile, 'TextType', 'string');
        for gi = 1:numel(gaits)
            r = T.gait == gaits{gi};
            if any(r), med(gaits{gi}) = T.medianFroude(find(r, 1)); end
        end
    end
    for gi = 1:numel(gaits)
        if ~isKey(med, gaits{gi})
            error(['GaitSel_ExportDryadMeanCycles: no speed threshold for %s. Run the R ' ...
                   'pipeline first so %s exists; recomputing it here would not match the ' ...
                   'figures.'], gaits{gi}, thrFile);
        end
    end

    binCell = {}; longRows = {};
    for gi = 1:numel(gaits)
        for sb = {'slow', 'fast'}
            sel = strcmp({S.gait}, gaits{gi});
            fr = [S.uBin] .^ 2;                       % Froude, as R compares it
            if strcmp(sb{1}, 'slow'), sel = sel & (fr <  med(gaits{gi}));
            else,                     sel = sel & (fr >= med(gaits{gi})); end
            grp = S(sel); if isempty(grp), continue; end
            us = [grp.u]; durs = [grp.dur];
            pct = linspace(0, 100, NP);
            % mean cycle duration for this bin maps the normalized cycle to a real-time base:
            % time_s = pct/100 * dur_mean_s. Individual steps vary (dur_sd_s), so this is the
            % mean-cycle timeline, not any one step.
            durMean = mean(durs); durSD = std(durs); tvec = pct / 100 * durMean;
            bin = struct('gait', gaits{gi}, 'speedbin', sb{1}, 'n', numel(grp), ...
                'u_median', median(us), 'u_min', min(us), 'u_max', max(us), ...
                'dur_mean_s', durMean, 'dur_sd_s', durSD, 'pct', pct, 'time_s', tvec);
            for ci = 1:numel(chans)
                Mx = cell2mat(arrayfun(@(s) s.cyc.(chans{ci}), grp, 'uni', 0)');
                mn = mean(Mx, 1, 'omitnan'); sd = std(Mx, 0, 1, 'omitnan');
                bin.([chans{ci} '_mean']) = mn; bin.([chans{ci} '_sd']) = sd;
                longRows{end+1} = table(repmat(string(gaits{gi}), NP, 1), repmat(string(sb{1}), NP, 1), ...
                    repmat(numel(grp), NP, 1), repmat(durMean, NP, 1), pct', tvec', ...
                    repmat(string(chans{ci}), NP, 1), mn', sd', ...
                    'VariableNames', {'gait','speedbin','n','dur_mean_s','pct','time_s','channel','mean','sd'}); %#ok<AGROW>
            end
            binCell{end+1} = bin; %#ok<AGROW>
        end
    end
    bins = [binCell{:}];
    save(fullfile(outDir, 'mean_cycles.mat'), 'bins', '-v7.3');
    writetable(vertcat(longRows{:}), fullfile(outDir, 'mean_cycles.csv'));
    fprintf('GaitSel_ExportDryadMeanCycles: %d gait x speed bins, %d channels -> %s\n', ...
        numel(bins), numel(chans), outDir);
end

function v = getOpt(o, n, d)
    if isstruct(o) && isfield(o, n) && ~isempty(o.(n)), v = o.(n); else, v = d; end
end
