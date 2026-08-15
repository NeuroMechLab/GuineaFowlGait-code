function [stepIdx, out] = detectStepsAccelMin(accelFA, fs, opts)
%DETECTSTEPSACCELMIN  Step detection from fore-aft CoM-acceleration braking minima.
%
%   [stepIdx, out] = detectStepsAccelMin(accelFA, fs)
%   [stepIdx, out] = detectStepsAccelMin(accelFA, fs, opts)
%
%   Detects one event per cycle at the MINIMUM of the fore-aft center-of-mass
%   acceleration - the instant of peak braking/deceleration, which occurs shortly
%   after each foot-strike. Because it keys on the CoM deceleration rather than on
%   the vertical force, it is robust across gaits (walking, running, hopping) and
%   works even when the force baseline is noisy or two feet overlap on one plate.
%
%   IMPORTANT CAVEAT: the braking minimum is a MID-CYCLE landmark, it does NOT
%   align with the start of stance (touchdown). Successive minima are a consistent
%   per-cycle boundary for cutting steps, but they do not give true stance-onset
%   timing. In this pipeline it is used only as the fallback cut source when a bout
%   has too few foot-marker touchdowns (see GaitSel_DetectGaitEvents).
%
%   INPUTS
%     accelFA : [n x 1] fore-aft CoM acceleration (the caller selects the fore-aft
%               column; the axis convention is yours). NaN-free.
%     fs      : sample rate (Hz).
%     opts    : struct with optional fields
%        .minPeakProminence  0.4  findpeaks MinPeakProminence, in the units of
%                                 accelFA (m/s^2). The main selectivity knob -
%                                 raise it to avoid extra peaks between steps.
%        .minCycleSec        []   if set, also imposes a minimum cycle time (s)
%                                 as findpeaks MinPeakDistance.
%
%   OUTPUTS
%     stepIdx : [nSteps x 1] sample indices of the braking minima, ascending.
%     out     : struct with .method, .eventType ('brakingMin'), .prominence,
%               .opts.
%
%   Port of the local CutStepsAcc (stepMode 1) from the GF turning pipeline
%   (ForcePlateData_StepDetectPerStepMeasures.m), as a pure function with a NaN
%   guard added. Requires the Signal Processing Toolbox (findpeaks). CC0 1.0.

if nargin < 3, opts = struct(); end
accelFA = accelFA(:);
if any(~isfinite(accelFA))
    error('detectStepsAccelMin:nonFinite', ...
          'accelFA must be finite (trim NaN padding to a clean window first).');
end

minProm     = getOpt(opts, 'minPeakProminence', 0.4);
minCycleSec = getOpt(opts, 'minCycleSec', []);

args = {'MinPeakProminence', minProm};
if ~isempty(minCycleSec)
    args = [args, {'MinPeakDistance', max(1, round(minCycleSec * fs))}];
end

% Minima of accelFA = peaks of the negated signal = peak fore-aft deceleration.
[~, pkLoc, ~, prom] = findpeaks(-accelFA, args{:});

stepIdx = pkLoc(:);
out = struct('method', 'accelMin', 'eventType', 'brakingMin', ...
             'prominence', prom(:), ...
             'opts', struct('minPeakProminence', minProm, 'minCycleSec', minCycleSec));
end


function v = getOpt(opts, name, default)
if isstruct(opts) && isfield(opts, name) && ~isempty(opts.(name))
    v = opts.(name);
else
    v = default;
end
end
