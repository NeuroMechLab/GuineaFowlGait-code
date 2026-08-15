function E = GaitSel_DetectGaitEvents(B, opts)
%GAITSEL_DETECTGAITEVENTS  Step segmentation + stance/flight timing for a bout.
%
%   E = GaitSel_DetectGaitEvents(B) segments the steps of a reconstructed bout B
%   and derives its stance / flight structure. Steps are cut at foot-marker
%   touchdowns (GaitSel_FootContacts), so each cycle starts at foot contact
%   (0 % = touchdown): running shows a single centred stance peak and walking its
%   double hump. When fewer than three foot contacts are found, the cut falls back
%   to the fore-aft CoM-acceleration braking minima (detectStepsAccelMin), a
%   gait-agnostic force/kinematics landmark. E.stepMethod records which was used.
%
%   The grounded/aerial structure comes from the summed vertical force: flight is
%   where it falls below a small fraction of body weight. Per-limb duty factor and
%   contact periods come from the foot-marker contacts (GaitSel_FootContacts).
%
%   E fields: stepIdx (cut samples used), stepMethod, stepPeriod, flightTime,
%   groundedFrac, hasFlight, stanceOn, stanceOff, dutyLimb (plate estimate),
%   dutyMarker (per-limb from markers), contactTimeMarker, footC, accelIdx, threshN.
%
%   opts: .flightThreshFrac (0.10 of BW), .minPeakProminence and .minStepMs (90)
%   for the accel-min fallback; foot-contact opts are forwarded to GaitSel_FootContacts.
%
%   See also: GaitSel_FootContacts, detectStepsAccelMin, GaitSel_PerStepStrideMeasures.

    if nargin < 2, opts = struct(); end
    g = 9.81; BW = B.mass*g; fHz = B.fHz; t = B.time;
    thr = getOpt(opts,'flightThreshFrac',0.10)*BW;

    % --- foot-marker contacts (per-limb; also the primary cycle-cut source) ---
    C = GaitSel_FootContacts(B, opts);

    % --- accel-minimum landmarks (fallback when < 3 foot contacts) -------
    aFA = filterForce(B.comAcc(:,2), fHz, struct('CutoffHz',20));
    prom = getOpt(opts,'minPeakProminence', 0.5*std(aFA,'omitnan'));
    minStep = getOpt(opts,'minStepMs',90)/1000;
    [accelIdx, sd] = detectStepsAccelMin(aFA, fHz, ...
        struct('minPeakProminence',prom, 'minCycleSec',minStep));

    % --- choose the cycle-cut points -------------------------------------
    footTD = sort(C.td(:));
    if ~isempty(footTD)                       % strictly increasing, no degenerate
        keep = [true; diff(footTD) >= round(0.03*fHz)];   % touchdowns >= 30 ms apart
        footTD = footTD(keep);
    end
    useFoot = numel(footTD) >= 3;
    if useFoot, stepIdx = footTD; else, stepIdx = accelIdx(:); end
    E.stepIdx    = stepIdx;
    E.stepMethod = string(char(iff(useFoot,'footTD','accelMin')));
    E.accelIdx   = accelIdx(:);
    E.footC      = C;
    E.brakeProm  = prom; E.accelOut = sd;

    % --- grounded / flight from summed vertical force --------------------
    Fz = filterForce(B.force(:,3), fHz, struct('CutoffHz',30));
    grounded = Fz > thr;

    nS = numel(stepIdx) - 1;
    [E.stepPeriod, E.flightTime, E.groundedFrac, E.stanceOn, E.stanceOff, E.dutyLimb] = ...
        deal(nan(max(nS,0),1));
    E.hasFlight = false(max(nS,0),1);
    for k = 1:nS
        win = stepIdx(k):stepIdx(k+1);
        gwin = grounded(win);
        E.stepPeriod(k)   = t(stepIdx(k+1)) - t(stepIdx(k));
        E.groundedFrac(k) = mean(gwin);
        E.flightTime(k)   = sum(~gwin)/fHz;
        E.hasFlight(k)    = any(~gwin) && E.flightTime(k) > 0.010;
        on = find(gwin,1,'first'); off = find(gwin,1,'last');
        if ~isempty(on); E.stanceOn(k) = win(on); E.stanceOff(k) = win(off); end
    end

    % --- supplementary per-limb duty factor from clean plate contacts ----
    E.dutyLimb = perLimbDuty(B, stepIdx, thr, BW, fHz);

    % --- per-limb duty factor + contact period from foot markers ---------
    E.contactTimeMarker = nan(max(nS,0),1);
    E.dutyMarker        = nan(max(nS,0),1);
    for k = 1:nS
        if k <= nS-1; sw = stepIdx(k):stepIdx(k+2); else; sw = stepIdx(max(k-1,1)):stepIdx(k+1); end
        dR = mean(C.downR(sw)); dL = mean(C.downL(sw));
        E.dutyMarker(k) = (dR + dL)/2;
        inStep = C.midIdx >= stepIdx(k) & C.midIdx < stepIdx(k+1);
        if any(inStep), E.contactTimeMarker(k) = max(C.stanceDur(inStep)); end
    end

    E.threshN = thr;
end

% -------------------------------------------------------------------- helpers
function dl = perLimbDuty(B, stepIdx, thr, BW, fHz)
    nS = numel(stepIdx)-1; dl = nan(max(nS,0),1);
    if ~isfield(B,'FzPlate'); return; end
    maxStance = 0.6; contacts = [];
    for p = 1:size(B.FzPlate,2)
        pk = max(B.FzPlate(:,p)); thp = max(thr, 0.20*pk);
        m = B.FzPlate(:,p) > thp;
        d = diff([0; m; 0]); s = find(d==1); e = find(d==-1)-1;
        for j = 1:numel(s)
            dur = (e(j)-s(j))/fHz;
            if dur >= 0.02 && dur <= maxStance, contacts = [contacts; s(j) e(j)]; end %#ok<AGROW>
        end
    end
    if isempty(contacts); return; end
    mid = mean(contacts,2);
    for k = 1:nS
        if k <= nS-1; strideP = (stepIdx(k+2)-stepIdx(k))/fHz; else; strideP = 2*(stepIdx(k+1)-stepIdx(k))/fHz; end
        inStep = mid >= stepIdx(k) & mid < stepIdx(k+1);
        if any(inStep), dl(k) = max((contacts(inStep,2)-contacts(inStep,1))/fHz) / strideP; end
    end
end

function out = iff(cond, a, b)
    if cond, out = a; else, out = b; end
end

function v = getOpt(o,n,d)
    if isstruct(o)&&isfield(o,n)&&~isempty(o.(n)); v=o.(n); else; v=d; end
end
