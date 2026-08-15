function GaitSelMulti_QCGate()
%GAITSELMULTI_QCGATE  QC-retention gate + steadiness labels (single source of truth).
%   Computes, from the tidy per-step/stride tables, the per-observation QC-pass flag
%   and the per-step energy-based steadiness class, and writes:
%       _RAnalysis/data/step_qcpass.csv     (boutID, stepIndex,   qcPass 0/1)
%       _RAnalysis/data/stride_qcpass.csv   (boutID, strideIndex, qcPass 0/1)
%       _RAnalysis/data/step_steadiness.csv (boutID, stepIndex,   accClass)  QC-passed steps
%
%   These are the SINGLE SOURCE OF TRUTH for the analyzed set: GaitSelMulti_ExportDryad
%   embeds qcPass/steadiness in the per-trial files, and the R pipeline (02_clean.R) reads
%   these flags rather than recomputing the gate, so the Dryad package, the mean traces,
%   and every R figure share one definition.
%
%   qcPass gate (matches 02_clean.R): finite positive Froude and meanSpeed; per-observation
%   work-energy residual WE_relresid <= 0.20 (missing passes); vertical CoM drift
%   driftRMS_vert_mm <= max(45 mm, DRIFT_RATE * trialDur_s) with DRIFT_RATE self-calibrated to
%   45 mm at the median-duration trial (missing passes); and not a freq/length/doubling outlier
%   (missing passes). trialDur_s is the per-bout sum of stepPeriod.
%   Steadiness (PRIMARY): net CoM mechanical energy change per unit distance as a fraction of
%   body weight, |dE_CoM/(m g L)| <= STEADY_GRADE (accelerating > +STEADY_GRADE, decelerating
%   < -STEADY_GRADE, else steady), an effective fore-aft grade. This is speed- and size-neutral,
%   because both the numerator and the distance scale grow with speed. Computed per stride and
%   inherited by its two steps; a step with no QC-passed parent stride falls back to its own
%   grade. Must match 02_clean.R.
%
%   See also: GaitSelMulti_ExportTidyCSV, GaitSelMulti_ExportDryad.

    WE_RELRESID_MAX = 0.20; DRIFT_MM = 45; STEADY_GRADE = 0.05; G = 9.81;
    P = projectPaths();
    St = readtable(fullfile(P.rData,'perStep_long.csv'),   'TextType','string');
    Sd = readtable(fullfile(P.rData,'perStride_long.csv'), 'TextType','string');

    % per-bout on-plate duration (sum of stepPeriod over finite positive steps) and the
    % self-calibrated drift rate (= DRIFT_MM at the median-duration bout)
    v = isfinite(St.stepPeriod) & St.stepPeriod > 0;
    [gb, bname] = findgroups(St.boutID(v));
    durByBout = splitapply(@sum, St.stepPeriod(v), gb);
    durMap = containers.Map(cellstr(bname), num2cell(durByBout));
    DRIFT_RATE = DRIFT_MM / median(durByBout);

    qStep   = qcFlag(St, durMap, DRIFT_MM, DRIFT_RATE, WE_RELRESID_MAX);
    qStride = qcFlag(Sd, durMap, DRIFT_MM, DRIFT_RATE, WE_RELRESID_MAX);

    if ~exist(P.rData,'dir'), mkdir(P.rData); end
    writetable(table(St.boutID, St.stepIndex,   qStep,   'VariableNames',{'boutID','stepIndex','qcPass'}), ...
               fullfile(P.rData,'step_qcpass.csv'));
    writetable(table(Sd.boutID, Sd.strideIndex, qStride, 'VariableNames',{'boutID','strideIndex','qcPass'}), ...
               fullfile(P.rData,'stride_qcpass.csv'));

    % ---- steadiness on the QC-passed set (energy grade = dE_CoM/(m g L)) ----
    cS = St(qStep==1, :); cD = Sd(qStride==1, :);
    strideEG    = cD.dE_CoM ./ (cD.mass_kg * G .* cD.strideLength);
    strideClass = classAcc(strideEG, STEADY_GRADE);
    strideKey   = cD.boutID + "|" + string(cD.strideIndex);
    sMap = containers.Map('KeyType','char','ValueType','char');
    for i = 1:height(cD), sMap(char(strideKey(i))) = char(strideClass(i)); end

    s0 = 2*floor((double(cS.stepIndex)-1)/2) + 1;              % parent stride index
    stepEG  = cS.dE_CoM ./ (cS.mass_kg * G .* cS.stepLength);  % fallback: step's own grade
    stepOwn = classAcc(stepEG, STEADY_GRADE);
    acc = strings(height(cS),1);
    for i = 1:height(cS)
        k = char(cS.boutID(i) + "|" + string(s0(i)));
        if isKey(sMap,k) && ~isempty(sMap(k)), acc(i) = string(sMap(k)); else, acc(i) = stepOwn(i); end
    end
    writetable(table(cS.boutID, cS.stepIndex, acc, 'VariableNames',{'boutID','stepIndex','accClass'}), ...
               fullfile(P.rData,'step_steadiness.csv'));

    fprintf('GaitSelMulti_QCGate: %d/%d steps, %d/%d strides pass QC; %d steady steps (%.0f%%).\n', ...
            sum(qStep), numel(qStep), sum(qStride), numel(qStride), sum(acc=="steady"), 100*mean(acc=="steady"));
end

% -------------------------------------------------------------------- helpers
function q = qcFlag(T, durMap, DRIFT_MM, DRIFT_RATE, WEMAX)
    n = height(T); td = zeros(n,1);
    for i = 1:n, k = char(T.boutID(i)); if isKey(durMap,k), td(i) = durMap(k); end, end  % missing -> 0
    we = T.WE_relresid; dr = T.driftRMS_vert_mm; ol = T.outlier;
    bound = max(DRIFT_MM, DRIFT_RATE .* td);
    pass = isfinite(T.Froude) & T.Froude > 0 & T.meanSpeed > 0 & ...
           (isnan(we) | we <= WEMAX) & ...
           (isnan(dr) | dr <= bound) & ...
           (isnan(ol) | ol == 0);
    q = double(pass);
end

% accelerating (x > thr) / decelerating (x < -thr) / steady (|x| <= thr); NaN -> "" (missing)
function c = classAcc(x, thr)
    c = strings(size(x));
    c(x >  thr) = "accelerating";
    c(x < -thr) = "decelerating";
    c(x >= -thr & x <= thr) = "steady";
end
