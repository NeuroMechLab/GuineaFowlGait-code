function GaitSelMulti_ExportCycleTraces(opts)
%GAITSELMULTI_EXPORTCYCLETRACES  Per-step cycle traces for the R hodograph figures.
%   Writes data/cycleTracesStep.csv: one STEP cycle per row-block, resampled to 100
%   points (0 % = touchdown), all channels dimensionless. This is the per-step,
%   unaveraged input the R pipeline registers and averages itself for Fig 4
%   (rotation sense), Fig 5 (steady hodographs), Fig 8 (transition hodographs) and
%   the transition-example supplement.
%
%   Built from the compiled Dryad per-trial series (GaitSel_DryadPackage_AllGF/
%   per_trial_timeseries), the same source as the mean-trace and mean-cycle exports,
%   so the reconstruction, fore-aft axis, step cuts and per-step indexing are
%   identical to the analyzed tables. Each step's samples startSample:endSample are
%   cut and resampled; channels are vertical and fore-aft force per body weight
%   (the series are already fore-aft-positive in travel), CoM kinetic and potential
%   energy fluctuation (about the cycle mean) per m*g*L0, and CoM fore-aft and
%   vertical velocity per sqrt(g*L0).
%
%   Each row carries the parent stride key so R can pair the two steps of a stride:
%   strideIndex = 2*floor((stepIndex-1)/2)+1 and stepInStride in {1,2}, the same
%   tiling the R code uses (stepIndex = strideIndex + stepInStride - 1).
    if nargin < 1, opts = struct(); end
    g = 9.81; nPts = 100; P = projectPaths(); addpath(P.helpers);
    root  = getOpt(opts, 'root', fullfile(P.root, 'GaitSel_DryadPackage_AllGF'));
    tsDir = fullfile(root, 'per_trial_timeseries');
    files = dir(fullfile(tsDir, '*.mat'));
    pct = (0:nPts-1)/(nPts-1)*100;

    stepRows = {};
    for f = 1:numel(files)
        T = load(fullfile(tsDir, files(f).name), 'trial').trial;
        st = T.steps; if isempty(st), continue; end
        m = T.mass_kg; L = T.L0_m; if ~(L > 0), continue; end
        BW = m*g; Vn = sqrt(g*L); n = numel(T.time_s);
        Bc = struct('com', T.com_m, 'comVel', T.comVel_ms, 'force', T.GRF_N, 'mass', m);
        for k = 1:height(st)
            a = st.startSample(k); b = min(st.endSample(k), n); if b - a < 6, continue; end
            win = a:b;
            cyc = cutCycle(Bc, win, BW, L, g, Vn, nPts);
            if isempty(cyc), continue; end
            si = double(st.stepIndex(k));
            s0 = 2*floor((si-1)/2) + 1;          % parent stride (odd step index)
            half = mod(si-1, 2) + 1;             % 1 or 2
            dur = T.time_s(b) - T.time_s(a);
            for p = 1:nPts
                stepRows{end+1} = struct('boutID', char(st.boutID(k)), 'strideIndex', s0, ...
                    'stepInStride', half, 'pct', pct(p), 'cycleDur_s', dur, ...
                    'Fz_BW', cyc.Fz_BW(p), 'Ffa_BW', cyc.Ffa_BW(p), 'KE_n', cyc.KE_n(p), ...
                    'PE_n', cyc.PE_n(p), 'vfa_n', cyc.vfa_n(p), 'vvert_n', cyc.vvert_n(p)); %#ok<AGROW>
            end
        end
    end
    if ~exist(P.rData,'dir'), mkdir(P.rData); end
    writetable(struct2table([stepRows{:}]), fullfile(P.rData,'cycleTracesStep.csv'));
    fprintf('GaitSelMulti_ExportCycleTraces: %d step-cycles -> cycleTracesStep.csv\n', numel(stepRows)/nPts);
end

% -------------------------------------------------------------------- helpers
% Cut one step window and resample every channel to nPts points on a
% percent-of-cycle axis. travelDir = +1 because the Dryad series are already
% oriented fore-aft-positive in the direction of travel.
function c = cutCycle(B, win, BW, L, g, Vn, nPts)
    c = [];
    if numel(win) < 6, return; end
    KE = 0.5*B.mass*sum(B.comVel(win,:).^2, 2); PE = B.mass*g*B.com(win,3);
    raw = struct('Fz_BW', B.force(win,3)/BW, 'Ffa_BW', B.force(win,2)/BW, ...
                 'KE_n', (KE-mean(KE))/(B.mass*g*L), 'PE_n', (PE-mean(PE))/(B.mass*g*L), ...
                 'vfa_n', B.comVel(win,2)/Vn, 'vvert_n', B.comVel(win,3)/Vn);
    fn = fieldnames(raw); x0 = linspace(1, numel(win), nPts);
    c = struct();
    for i = 1:numel(fn), c.(fn{i}) = interp1(1:numel(win), raw.(fn{i}), x0)'; end
end

function v = getOpt(o,n,d), if isstruct(o)&&isfield(o,n)&&~isempty(o.(n)), v=o.(n); else, v=d; end, end
