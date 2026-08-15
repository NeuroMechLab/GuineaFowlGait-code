function C = GaitSel_FootContacts(B, opts)
%GAITSEL_FOOTCONTACTS  Per-limb foot contact periods from the foot markers.
%
%   C = GaitSel_FootContacts(B) detects the ground-contact intervals of each foot
%   (right, left) from the foot markers of an imported bout B, independent of the
%   plate-force thresholding that merges walking footfalls. The foot-height signal is
%   first cleaned (physically implausible marker spikes and robust MAD outliers removed,
%   gaps filled and lightly low-pass smoothed), because untracked-marker extrapolation
%   spikes otherwise corrupt the contact baseline and cause detection to fail on slow
%   walking. A foot is then in contact when it is low (near the ground): contact = cleaned
%   foot height below an ADAPTIVE baseline (5th percentile of the cleaned signal) plus a
%   fixed tolerance (default 3 cm). Short dropouts within a contact are bridged and runs
%   shorter than the minimum contact are dropped. This reproduces the documented walking
%   duty factor (about 0.66) and, for walking, makes the union of the two feet's contacts
%   equal 1 (no flight), consistent with the force. The markers thus split the ground-contact
%   time between limbs, including the double support of walking that the summed five-plate
%   force cannot resolve, and the per-foot touchdowns (C.td) are the cycle-cut source used by
%   GaitSel_DetectGaitEvents. The union of the contacts is reported against the force grounded
%   fraction as a consistency check (they should agree).
%
%   The foot markers are on the force time base (GaitSelMulti_ImportBout resamples
%   them), so contacts align with the reconstructed CoM and events.
%
%   opts: .zTolM (contact height band, default 0.03 m), .minContactMs (default 40),
%         .bridgeMs (bridge in-contact dropouts up to this long, default 20).
%
%   OUTPUT struct C: downR, downL (logical masks over the force time base),
%   td, to, side, stanceDur, midIdx (per contact), zTol, unionFrac (marker),
%   groundedTarget (force grounded fraction).
%
%   See also: GaitSel_DetectGaitEvents, GaitSel_PerStepStrideMeasures.

    if nargin < 2, opts = struct(); end
    fHz = B.fHz;
    zTol = getOpt(opts,'zTolM', 0.03);
    minC = round(getOpt(opts,'minContactMs',40)/1000 * fHz);
    bridgeN = round(getOpt(opts,'bridgeMs',20)/1000 * fHz);   % close short in-contact dropouts

    % The foot tracks arrive already spike-cleaned from GaitSelMulti_ImportBout, so they are
    % thresholded as they are. Cleaning again here would low-pass the height a second time and
    % put the contact detector on a different signal from the virtual leg and the stored series.
    zR = B.footR(:,3);
    zL = B.footL(:,3);
    gR = quantile(zR,0.05); gL = quantile(zL,0.05);          % adaptive baseline on cleaned signal

    % force-derived grounded fraction (consistency reference)
    thr = 0.10*B.mass*9.81;
    grounded = filterForce(B.force(:,3), fHz, struct('CutoffHz',30)) > thr;
    target = mean(grounded);

    downR = zR < gR + zTol;  downL = zL < gL + zTol;
    unionFrac = mean(downR | downL);

    % contact intervals per foot for QC / stance durations
    td=[]; to=[]; side=strings(0); dur=[]; mid=[];
    for s = ["R","L"]
        dn = (s=="R").*downR + (s=="L").*downL > 0;
        runs = contiguousRuns(dn, minC, bridgeN);
        for r = 1:size(runs,1)
            td(end+1,1)=runs(r,1); to(end+1,1)=runs(r,2); %#ok<AGROW>
            side(end+1,1)=s; dur(end+1,1)=(runs(r,2)-runs(r,1))/fHz; %#ok<AGROW>
            mid(end+1,1)=round((runs(r,1)+runs(r,2))/2); %#ok<AGROW>
        end
    end
    [td,ord]=sort(td); to=to(ord); side=side(ord); dur=dur(ord); mid=mid(ord);
    C = struct('downR',downR,'downL',downL,'td',td,'to',to,'side',side, ...
               'stanceDur',dur,'midIdx',mid,'zTol',zTol,'unionFrac',unionFrac, ...
               'groundedTarget',target);
end

function runs = contiguousRuns(b, minLen, bridgeN)
    b = b(:) > 0;
    if nargin >= 3 && bridgeN > 0          % bridge short false-gaps within a contact
        nb = ~b; d = diff([0; nb; 0]);     % false-runs of b, in b-coordinates
        gs = find(d==1); ge = find(d==-1)-1;
        for k = 1:numel(gs)
            if (ge(k)-gs(k)+1) <= bridgeN, b(gs(k):ge(k)) = true; end
        end
    end
    d = diff([0; b; 0]);
    s = find(d==1); e = find(d==-1)-1;
    keep = (e-s+1) >= minLen;
    runs = [s(keep) e(keep)];
end

% Foot-track cleaning happens once, in GaitSelMulti_ImportBout, via cleanFootTrack. Every
% consumer therefore sees identical feet: this detector, the virtual-leg geometry in
% GaitSel_PerStepStrideMeasures, the CoM-offset fit, and the published per-trial series.

function v = getOpt(o,n,d)
    if isstruct(o)&&isfield(o,n)&&~isempty(o.(n)); v=o.(n); else; v=d; end
end
