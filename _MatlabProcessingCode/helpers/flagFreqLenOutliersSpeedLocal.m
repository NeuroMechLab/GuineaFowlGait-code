function flag = flagFreqLenOutliersSpeedLocal(freq, len, speed, k)
%FLAGFREQLENOUTLIERSSPEEDLOCAL  Speed-conditional step frequency/length outlier flag.
%   flag = flagFreqLenOutliersSpeedLocal(freq, len, speed, k) flags steps whose frequency OR
%   length is anomalous RELATIVE TO OTHER STEPS AT THE SAME SPEED, rather than against a global
%   median. Because forward speed = step length x step frequency, and the dataset spans walking
%   through fast running, a global log-MAD filter (logMadFlag) has a running-dominated median and
%   clips the legitimately low-frequency, long-step slow-walk tail as "outliers". Judging each step
%   within its dimensionless-speed neighbourhood removes that bias while still catching genuine
%   mis-cut cycles (which are anomalous even among their speed peers).
%
%   Speed is split into ~12 quantile bins; within each populated bin (>= 8 steps) a log-MAD flag
%   (threshold k) is applied to frequency and to length. Sparse bins and degenerate inputs fall
%   back to the global logMadFlag so the filter never fails open. Returns a logical column vector.
%   See GaitSelMulti_ExportTidyCSV, flagSpeedDoubling, logMadFlag.
    freq = freq(:); len = len(:); speed = speed(:);
    n = numel(freq); flag = false(n,1);
    ok = isfinite(freq) & freq > 0 & isfinite(len) & len > 0 & isfinite(speed) & speed > 0;
    gf = logMadFlag(freq, k); gl = logMadFlag(len, k);          % global fallback
    if nnz(ok) < 24
        flag = gf | gl; return
    end
    edges = unique(quantile(speed(ok), linspace(0,1,13)));
    if numel(edges) < 3
        flag = gf | gl; return
    end
    [~,~,bin] = histcounts(speed, edges);
    lf = log(freq); ll = log(len);
    flag(~ok) = true;                                          % non-finite/non-positive: flag
    for b = 1:max(bin)
        sel = ok & bin == b;
        if nnz(sel) < 8
            flag(sel) = gf(sel) | gl(sel);                     % too few peers: use global
            continue
        end
        mf = median(lf(sel)); sf = 1.4826*median(abs(lf(sel)-mf));
        ml = median(ll(sel)); sl = 1.4826*median(abs(ll(sel)-ml));
        flag(sel) = abs(lf(sel)-mf) > k*max(sf,eps) | abs(ll(sel)-ml) > k*max(sl,eps);
    end
end
