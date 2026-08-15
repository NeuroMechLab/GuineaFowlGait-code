function z = cleanFootSignal(z, fHz)
%CLEANFOOTSIGNAL  Remove marker spikes, fill gaps, and lightly smooth a foot-height channel.
%   Untracked foot markers are filled/extrapolated to physically impossible heights at the bout
%   edges; those spikes corrupt contact-baseline percentiles and foot-relative measures. Set
%   implausible samples (|height| > 0.4 m) and robust MAD outliers (> 6 robust s.d. from the
%   median) to missing, fill by shape-preserving interpolation holding the ends, then apply a
%   zero-phase 20 Hz low-pass (foot motion is well below this). Used by GaitSel_FootContacts
%   The guinea fowl pipeline does not call this: foot tracks are cleaned once at import by
%   cleanFootTrack, which applies the same criterion across all three axes together. This
%   single-channel form is kept for readers that carry only a height signal. See also
%   cleanFootTrack.
    z = z(:); z(abs(z) > 0.4) = NaN;
    ok = isfinite(z);
    if nnz(ok) < 10, z = fillmissing(z, 'constant', median(z,'omitnan')); return; end
    m = median(z(ok)); s = 1.4826*median(abs(z(ok)-m));
    z(abs(z-m) > 6*max(s,eps)) = NaN;
    z = fillmissing(z, 'pchip', 'EndValues', 'nearest');
    z = filterForce(z, fHz, struct('CutoffHz', 20));
end
