function flag = flagSpeedDoubling(freq, speed, fracThresh)
%FLAGSPEEDDOUBLING  Speed-conditional mis-cut-cycle flag (SYMMETRIC).
%   flag = flagSpeedDoubling(freq, speed, fracThresh) flags step/stride cycles whose
%   frequency is inconsistent with the local norm at their forward speed. At a given
%   dimensionless speed a cycle's frequency should sit near the speed-local median. A
%   missed-touchdown MERGED cycle has ~half that frequency (same speed, ~2x period); a
%   spurious extra-touchdown SPLIT cycle has ~twice it (~half period). Both are detection
%   artifacts that are invisible to the speed-blind global log-MAD filter, because forward
%   speed = step length x step frequency is preserved. This flags EITHER tail: a
%   log-frequency more than |log(fracThresh)| from the speed-bin median (fracThresh = 0.65
%   flags below 0.65x or above 1.54x the local median). Speed is binned into ~12 quantile
%   groups; the median is robust to the (minority) mis-cut cycles within a bin. Returns a
%   logical column vector the size of freq; all-false if there are too few valid cycles
%   (< 24) or too few distinct speed bins to form a stable local median.
%
%   Used by GaitSelMulti_ExportTidyCSV (the outlier flag).
    freq = freq(:); speed = speed(:);
    flag = false(size(freq));
    ok = isfinite(freq) & freq > 0 & isfinite(speed) & speed > 0;
    if nnz(ok) < 24, return; end
    edges = unique(quantile(speed(ok), linspace(0,1,13)));
    if numel(edges) < 3, return; end
    [~,~,bin] = histcounts(speed, edges);
    lf = log(freq); tol = abs(log(fracThresh));
    for b = 1:max(bin)
        sel = ok & bin == b;
        if nnz(sel) < 8, continue; end
        flag(sel) = abs(lf(sel) - median(lf(sel))) > tol;
    end
end
