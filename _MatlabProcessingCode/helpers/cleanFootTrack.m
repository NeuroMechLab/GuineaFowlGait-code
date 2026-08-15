function foot = cleanFootTrack(foot, fHz)
%CLEANFOOTTRACK  Remove marker spikes from a 3-D foot-marker track, consistently across axes.
%   foot is n-by-3 [medio-lateral, fore-aft, vertical]. An untracked/extrapolated foot marker
%   spikes in ALL axes at the same samples (e.g. a physically impossible 5 m foot height at a
%   bout edge), which corrupts foot-relative measures such as the virtual leg length. Bad samples
%   are detected from the vertical channel (implausible height > 0.4 m, or a robust MAD outlier)
%   and set to missing in ALL THREE axes, then each axis is gap-filled (shape-preserving, ends
%   held) and lightly low-pass smoothed. This keeps the axes time-consistent (unlike cleaning each
%   axis independently).
%
%   This is called ONCE per bout, in GaitSelMulti_ImportBout, so the contact detector, the
%   virtual-leg geometry, the CoM-offset fit and the published per-trial series all measure the
%   same feet. Do not call it again downstream: a second pass would low-pass the track twice and
%   reintroduce the split it exists to remove. See also cleanFootSignal.
    z = foot(:,3);
    bad = abs(z) > 0.4;
    ok = isfinite(z) & ~bad;
    if nnz(ok) >= 10
        m = median(z(ok)); s = 1.4826*median(abs(z(ok)-m));
        bad = bad | (abs(z - m) > 6*max(s,eps));
    end
    for ax = 1:3
        c = foot(:,ax); c(bad) = NaN;
        if nnz(isfinite(c)) < 10
            foot(:,ax) = fillmissing(foot(:,ax), 'constant', median(foot(:,ax),'omitnan'));
            continue;
        end
        c = fillmissing(c, 'pchip', 'EndValues', 'nearest');
        foot(:,ax) = filterForce(c, fHz, struct('CutoffHz', 20));
    end
end
