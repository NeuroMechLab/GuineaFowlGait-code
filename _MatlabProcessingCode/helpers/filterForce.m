function y = filterForce(x, fs, opts)
%FILTERFORCE  Zero-phase Butterworth low-pass for ground-reaction-force signals.
%   y = filterForce(x, fs) low-passes each column of X (force / moment channels,
%   one signal per column) with a 2nd-order Butterworth filter applied forward and
%   backward via FILTFILT (zero phase; effective 4th order), at a 30 Hz cutoff.
%
%   y = filterForce(x, fs, opts) with opts.CutoffHz / opts.Order overrides.
%
%   USE THIS FOR FORCE/MOMENT CHANNELS ONLY. Force-plate signals do not drop out
%   and a value of exactly 0 N is genuine (no load), so X must be finite and is
%   filtered as-is with no gap handling.
%
%   Standard preset (cf. Winter, Biomechanics and Motor Control of Human
%   Movement):  CutoffHz = 30 Hz,  Order = 2.
%
%   Vendored from the E117 Shared_Code/filterForce.m (adapted to an opts struct).
%   Requires the Signal Processing Toolbox (butter, filtfilt). CC0 1.0.

if nargin < 3, opts = struct(); end
cutoffHz = getOpt(opts, 'CutoffHz', 30);
order    = getOpt(opts, 'Order', 2);

if ~isscalar(fs) || ~(fs > 0)
    error('filterForce:fs', 'fs must be a positive scalar (Hz).');
end
if cutoffHz >= fs/2
    error('filterForce:cutoffAboveNyquist', ...
        'CutoffHz (%g Hz) must be below the Nyquist frequency (%g Hz).', cutoffHz, fs/2);
end
if any(~isfinite(x(:)))
    error('filterForce:nonFinite', ...
        ['Input contains NaN/Inf. filterForce is for clean force channels; ' ...
         'trim or gap-fill before filtering.']);
end

[b, a] = butter(order, cutoffHz / (fs/2), 'low');
y = filtfilt(b, a, x);   % operates column-wise on matrices
end


function v = getOpt(opts, name, default)
if isstruct(opts) && isfield(opts, name) && ~isempty(opts.(name))
    v = opts.(name);
else
    v = default;
end
end
