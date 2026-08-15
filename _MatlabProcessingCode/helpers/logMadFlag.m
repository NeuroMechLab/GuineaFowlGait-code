function flag = logMadFlag(x, k)
%LOGMADFLAG  Robust log-space median-absolute-deviation outlier flag.
%   flag = logMadFlag(x, k) flags entries of x whose log deviates from the median log by more
%   than k robust standard deviations (1.4826*MAD). Non-finite or non-positive entries are
%   flagged (they cannot be a valid frequency/length). Used by the step frequency/length outlier
%   filters. NOTE: applied globally this penalizes the naturally low-frequency slow-walk tail of a
%   bimodal walk+run distribution; prefer flagFreqLenOutliersSpeedLocal for speed-covarying
%   quantities. See GaitSelMulti_ExportTidyCSV.
    x = x(:); ok = isfinite(x) & x > 0; lx = log(x);
    m = median(lx(ok)); s = 1.4826*median(abs(lx(ok)-m));
    flag = false(size(x));
    flag(ok) = lx(ok) < (m-k*s) | lx(ok) > (m+k*s);
    flag(~ok) = true;
end
