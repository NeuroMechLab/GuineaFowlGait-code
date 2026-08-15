function [S, T] = GaitSelMulti_ExportTidyCSV()
%GAITSELMULTI_EXPORTTIDYCSV  Normalise the per-step/stride measures, write the tidy tables.
%   Reads SavedResults/perStepMeasures_multi.mat, derives each (study-tagged) bird's L0
%   from its touchdown leg length over every non-aerial step (within-sample allometry as the
%   fallback for a bird with too few), adds the dimensionless columns, flags freq/length outliers, and writes the tidy tables
%   under both the _multi names and the canonical names the R pipeline reads (01_load.R):
%       _RAnalysis/data/perStride_long_multi.csv  ->  perStride_long.csv
%       _RAnalysis/data/perStep_long_multi.csv    ->  perStep_long.csv
%       _RAnalysis/data/morphology_multi.csv       ->  morphology.csv
%
%   Normalisers (dynamic similarity): force = m g, length = L0, velocity =
%   sqrt(g L0), work = m g L0, Froude = v^2/(g L0). Study/dataset/colourCode
%   columns pass through unchanged.
%
%   See also: GaitSelMulti_BatchProcess, GaitSelMulti_ExportDryad.

    P = projectPaths(); addpath(P.helpers); g = 9.81;
    load(fullfile(P.saved,'perStepMeasures_multi.mat'),'S','T','massBird');

    % ---- L0: median touchdown leg length over every NON-AERIAL step ------
    % The subset is every step with no aerial phase, walk and grounded running together, read from
    % hasFlight. hasFlight comes from the summed vertical force alone, so the subset carries no
    % dependence on the gait classification, which matters because L0 sets the dimensionless speed
    % u, and u is the abscissa of the rotation-sense crossover fit.
    %
    % The non-aerial subset also reaches the five-step minimum for more bird-sessions, since the
    % Blum cohorts were run fast and several of their birds produced no walk-class step.
    %
    % L0 is derived per BIRD-SESSION, on the same birdCode|YYYY-MM key as the mass and the CoM
    % offset. Touchdown leg length is measured from the CoM proxy to the toe, so it moves with
    % both the marker placement and the offset fitted for that session; pooling a bird's sessions
    % would mix two CoM definitions into one length scale.
    %
    % Where a session holds too few non-aerial steps the estimate falls back to the bird, then to
    % a within-sample allometry of the measured L0 on body mass, and the branch taken is recorded
    % per session in morphology.csv so no fallback is silent. The allometry is fitted on this
    % sample rather than taken as 0.20*m^(1/3): the geometric form is a different quantity, an
    % isometric reference length rather than a measured leg, and mixing the two would put parts of
    % the sample on different length scales.
    T.session = sessionKey(T.bird, T.dateCode);
    S.session = sessionKey(S.bird, S.dateCode);
    sessions = unique(T.session);
    L0 = containers.Map('KeyType','char','ValueType','double');
    Lnorm = containers.Map('KeyType','char','ValueType','double');
    L0src = containers.Map('KeyType','char','ValueType','char');
    okLen = isfinite(T.legLen_TD) & T.legLen_TD>0.10 & T.legLen_TD<0.40;
    nonAerial = T.hasFlight == 0;
    MIN_LEN_STEPS = 5;
    mHave = []; lHave = []; pend = {};
    for i = 1:numel(sessions)
        b = char(sessions(i)); mkg = massBird(b);
        Lnorm(b) = 0.20 * mkg^(1/3);
        sel = T.session==sessions(i) & nonAerial & okLen;
        if nnz(sel) >= MIN_LEN_STEPS
            L0(b) = median(T.legLen_TD(sel)); L0src(b) = 'session_nonAerial_median';
            mHave(end+1,1) = mkg; lHave(end+1,1) = L0(b); %#ok<AGROW>
        else
            L0(b) = NaN; L0src(b) = 'pending';
            pend{end+1} = b; %#ok<AGROW>
            fprintf('  %s: only %d non-aerial steps, L0 from the bird\n', b, nnz(sel));
        end
    end

    % Second choice: the same BIRD's non-aerial steps pooled over its sessions. A session run
    % fast can hold too few grounded steps to measure a leg from, and the leg measured in that
    % bird's other sessions is the same leg: for the birds with more than one measured session the
    % session estimates agree far more closely with each other than birds do with one another. The
    % cross-bird allometry below discards the individual, so it is the last resort rather than the
    % first fallback. morphology.csv records which branch set each value.
    if ~isempty(pend)
        still = {};
        for i = 1:numel(pend)
            b = pend{i};
            sel = T.bird == extractBefore(string(b),"|") & nonAerial & okLen;
            if nnz(sel) >= MIN_LEN_STEPS
                L0(b) = median(T.legLen_TD(sel)); L0src(b) = 'bird_nonAerial_median';
                fprintf('  %s: L0 %.4f m from the bird''s %d non-aerial steps across sessions\n', ...
                        b, L0(b), nnz(sel));
            else
                still{end+1} = b; %#ok<AGROW>
            end
        end
        pend = still;
    end

    if ~isempty(pend)
        if numel(mHave) >= 4
            p = polyfit(log(mHave), log(lHave), 1);
            for i = 1:numel(pend)
                b = pend{i};
                L0(b) = exp(p(2)) * massBird(b)^p(1);
                L0src(b) = 'withinSample_allometry';
                warning('GaitSelMulti:L0allometry', ...
                    '%s: L0 %.4f m from the within-sample allometry (exponent %.3f)', ...
                    b, L0(b), p(1));
            end
        else
            for i = 1:numel(pend)
                b = pend{i}; L0(b) = Lnorm(b); L0src(b) = 'geometric_lastResort';
                warning('GaitSelMulti:L0geom','%s: L0 fell back to 0.20*m^(1/3)', b);
            end
        end
    end

    S = addNorm(S, L0, massBird, g, 'stride');
    T = addNorm(T, L0, massBird, g, 'step');
    T = addLegStiffness(T, g);      % needs L0_m/BW_N, so it follows addNorm

    % Outlier flag = a SPEED-CONDITIONAL log-MAD frequency/length filter (k=3.0) UNION a
    % SYMMETRIC speed-conditional doubling filter. Both are judged within dimensionless-speed
    % neighbourhoods rather than globally, because forward speed = step length x step frequency
    % and the dataset spans walking through fast running: a GLOBAL log-MAD frequency filter has a
    % running-dominated median and clips the legitimately low-frequency, long-step slow-walk tail
    % as "outliers" (a speed-selection bias). flagFreqLenOutliersSpeedLocal judges each step
    % against its speed peers, so genuine mis-cuts are still caught while slow walks are not
    % penalised for being slow. The doubling filter catches merged (half-frequency) and split
    % (double-frequency) mis-cut cycles that keep the same forward speed: with DBL_FRAC = 0.65 it
    % flags any cycle whose log frequency departs from its speed-local median by more than
    % |log(0.65)| in EITHER direction.
    OUTLIER_K = 3.0; DBL_FRAC = 0.65;
    T.outlier = double(logical(flagFreqLenOutliersSpeedLocal(T.stepFreq,   T.stepLength,   T.meanSpeed_n, OUTLIER_K)) | ...
                       flagSpeedDoubling(T.stepFreq,   T.meanSpeed_n, DBL_FRAC));
    S.outlier = double(logical(flagFreqLenOutliersSpeedLocal(S.strideFreq, S.strideLength, S.meanSpeed_n, OUTLIER_K)) | ...
                       flagSpeedDoubling(S.strideFreq, S.meanSpeed_n, DBL_FRAC));
    fprintf('  outliers: %d/%d steps, %d/%d strides (speed-local log-MAD freq/len UNION speed-conditional doubling).\n', ...
            sum(T.outlier), height(T), sum(S.outlier), height(S));

    if ~exist(P.rData,'dir'); mkdir(P.rData); end
    writetable(S, fullfile(P.rData,'perStride_long_multi.csv'));
    writetable(T, fullfile(P.rData,'perStep_long_multi.csv'));

    nb = numel(sessions);
    M = table('Size',[nb 9], 'VariableTypes', ...
        {'string','string','string','string','double','double','double','double','double'}, ...
        'VariableNames', {'session','bird','study','sessionMonth','mass_kg','L0_m','Lnorm_m','nSteps','nStrides'});
    M.L0_source = strings(nb,1);   % which branch set L0, so a fallback is never silent
    for i = 1:nb
        b = char(sessions(i));
        j = find(T.session==sessions(i),1);
        M.session(i)=string(b); M.bird(i)=T.bird(j); M.study(i)=T.study(j);
        M.sessionMonth(i)=extractAfter(string(b),"|");
        M.mass_kg(i)=massBird(b);
        M.L0_m(i)=L0(b); M.Lnorm_m(i)=Lnorm(b); M.L0_source(i)=string(L0src(b));
        M.nSteps(i)=nnz(T.session==sessions(i)); M.nStrides(i)=nnz(S.session==sessions(i));
    end
    writetable(M, fullfile(P.rData,'morphology_multi.csv'));

    % Promote the multi tables to the canonical names the R pipeline reads (01_load.R reads
    % perStep_long.csv / perStride_long.csv / morphology.csv). The unified multi-study dataset IS
    % the analyzed dataset, so this makes it the canonical hand-off automatically (no manual copy).
    copyfile(fullfile(P.rData,'perStep_long_multi.csv'),   fullfile(P.rData,'perStep_long.csv'));
    copyfile(fullfile(P.rData,'perStride_long_multi.csv'), fullfile(P.rData,'perStride_long.csv'));
    copyfile(fullfile(P.rData,'morphology_multi.csv'),     fullfile(P.rData,'morphology.csv'));
    fprintf('GaitSelMulti_ExportTidyCSV: %d strides, %d steps, %d bird-sessions from %d birds (also promoted to canonical names).\n', ...
            height(S), height(T), nb, numel(unique(T.bird)));
end

% -------------------------------------------------------------------- helpers

function k = sessionKey(bird, dateCode)
%SESSIONKEY  birdCode|YYYY-MM, the key mass, L0 and the CoM offset all share.
    t = regexprep(string(dateCode), '\D', '');
    assert(all(strlength(t) >= 6), 'GaitSelMulti:badDateCode', ...
           'a dateCode has no YYYYMM to key the session on');
    k = string(bird) + "|" + extractBetween(t,1,4) + "-" + extractBetween(t,5,6);
end
function X = addNorm(X, L0, massBird, g, kind)
    b = arrayfun(@(x) char(x), X.session, 'uni', 0);
    L = cellfun(@(k) L0(k), b); m = cellfun(@(k) massBird(k), b);
    X.L0_m = L; X.BW_N = m*g; X.mass_kg = m;
    common = {'meanSpeed','velocity'; 'contactTime','time'; 'flightTime','time'; ...
              'accelForeAft','acceleration'; 'hodoArea',[0 2 -2]};
    if strcmp(kind,'stride')
        dimMap = [common; {
            'strideLength','length'; 'strideFreq','frequency'; 'stridePeriod','time'; ...
            'dKE','work'; 'dPE','work'; 'W_CoM_pos','work'; 'W_CoM_neg','work'; ...
            'W_CoM_net','work'; 'Wf_horiz','work'; 'Wv_vert','work'; 'Wtot_com','work'; ...
            'KE_range','work'; 'PE_range','work'; 'vertExcursion','length'}];
    else
        dimMap = [common; {
            'speed_TD','velocity'; 'stepLength','length'; 'stepPeriod','time'; 'stepFreq','frequency'; ...
            'legLen_TD','length'; 'impulseForeAft','impulse'; ...
            'W_CoM_net','work'; 'dE_CoM','work'; 'dKE_step','work'; 'dPE_step','work'; ...
            'W_CoM_pos','work'; 'W_CoM_neg','work'; 'KE_range','work'; 'PE_range','work'; ...
            'vertExcursion','length'}];
    end
    dimMap = dimMap(ismember(dimMap(:,1), X.Properties.VariableNames), :);
    X = scaleByDimension(X, dimMap, 'mass_kg', 'L0_m', 'g', g);
    X.Froude = X.meanSpeed_n.^2;
end

function X = addLegStiffness(X, g)
%ADDLEGSTIFFNESS  Dimensionless leg-spring stiffness per step (Blum et al. 2009 method C).
%
%   Adds two columns to the per-step tidy table:
%     legCompress_n  leg compression as a fraction of touchdown leg length, dL/L0
%     kLeg_n         dimensionless leg stiffness k L0/(m g), i.e. body weights per leg
%                    length (the Geyer, Seyfarth & Blickhan 2006 convention)
%
%   Leg stiffness is k = Fmax/dL with Fmax the peak vertical force. The leg compression dL
%   is recovered by inverting a sinusoidal vertical force profile over the contact period
%   (Blum, Lipfert & Seyfarth 2009, their method C):
%
%       dL = L0 + (Fmax/m)(tc/pi)^2 - (g/8) tc^2 - L0 sin(alpha_TD)
%
%   with tc the per-limb contact time and alpha_TD the leg angle from the horizontal at
%   touchdown, which is exactly the legAngle_TD convention used here. Blum et al. compared
%   five leg-compression conventions against bipedal spring-mass predictions: method C
%   agrees with the model where the symmetric-trajectory convention (McMahon & Cheng 1990;
%   Farley, Glasheen & McMahon 1993) underestimates stiffness. Their method E is not used
%   because it infers Fmax from duty factor and pins it at pi/2 body weights once duty
%   factor reaches 0.5, which every grounded step does.
%
%   Steps whose inverted compression is non-physical (dL <= 0, or more than 0.4 L0) get
%   NaN rather than a meaningless ratio. The surviving per-step distribution is
%   right-skewed, so summaries use median and interquartile range (R 22_table_gait4_summary).
%
%   Leg stiffness is a MODEL-BOUND descriptor: it exists only relative to the spring-mass
%   template it is derived from, and where two limbs share the net force in double support
%   it is an effective two-limb stiffness, not a single-limb property (McMahon 1985; Farley
%   et al. 1993; Geyer et al. 2006). No effective vertical stiffness is emitted: it
%   corresponds to no physical spring in the model, and its touchdown-referenced
%   displacement degenerates in a vaulting gait, where the CoM rises after touchdown.
%

    need = {'peakVertForce_BW','contactTime','legAngle_TD','L0_m','BW_N','mass_kg'};
    if ~all(ismember(need, X.Properties.VariableNames))
        X.legCompress_n = nan(height(X),1); X.kLeg_n = nan(height(X),1);
        warning('GaitSelMulti:legStiffness','missing inputs; kLeg_n set to NaN');
        return
    end
    L0    = X.L0_m;
    tc    = X.contactTime;
    Fmax  = X.peakVertForce_BW .* X.BW_N;                 % N
    alpha = X.legAngle_TD * pi/180;                        % from horizontal at touchdown
    dL    = L0 + (Fmax ./ X.mass_kg) .* (tc/pi).^2 - (g/8)*tc.^2 - L0 .* sin(alpha);
    frac  = dL ./ L0;
    ok    = isfinite(frac) & frac > 0 & frac < 0.40;
    X.legCompress_n = nan(height(X),1); X.legCompress_n(ok) = frac(ok);
    X.kLeg_n        = nan(height(X),1);
    X.kLeg_n(ok)    = X.peakVertForce_BW(ok) ./ frac(ok);
    fprintf('  leg stiffness (Blum 2009 method C): kLeg_n on %d/%d steps, median %.1f BW/L0.\n', ...
            nnz(ok), height(X), median(X.kLeg_n(ok),'omitnan'));
end

% The outlier kernels are standalone helpers (helpers/flagSpeedDoubling.m,
% helpers/flagFreqLenOutliersSpeedLocal.m, helpers/logMadFlag.m), on the path via
% projectPaths. The speed-local judging is explained in the outlier-flag comment above.
