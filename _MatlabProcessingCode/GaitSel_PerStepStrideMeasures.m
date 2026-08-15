function [stepRows, strideRows] = GaitSel_PerStepStrideMeasures(B, E)
%GAITSEL_PERSTEPSTRIDEMEASURES  Per-step and per-stride measures for one bout.
%
%   [stepRows, strideRows] = GaitSel_PerStepStrideMeasures(B, E) reduces a
%   reconstructed bout B (GaitSel_ReconstructCoM) and its gait events E
%   (GaitSel_DetectGaitEvents) to one row per step (forces, timing, per-step
%   speed and fore-aft acceleration, virtual leg) and one row per stride (CoM
%   energy exchange: kinetic and gravitational potential energy fluctuations,
%   pendular recovery, congruity, and an objective gait classification).
%
%   Steps span the cut samples E.stepIdx (foot-marker touchdowns); a stride is two
%   steps. Stance within a step is the grounded window (E.stanceOn:E.stanceOff). Energy
%   exchange follows Cavagna, Heglund & Taylor (1977): forward kinetic energy
%   Ekf = 1/2 m (vml^2 + vfa^2); vertical energy Ev = m g h + 1/2 m vv^2;
%   recovery R = (Wf + Wv - Wtot)/(Wf + Wv), each W the sum of positive
%   increments over the stride. Congruity is the fraction of the stride in which
%   Ekf and Ev change in the same direction (high = in-phase bouncing).
%
%   Gait is classified by two threshold-free mechanistic features: aerial phase
%   present (aerial run), else grounded and split by the hodograph rotation sense
%   (sign of the CoM velocity-loop signed area): pendular/vaulting = walk,
%   bouncing = grounded run. Recovery, congruity, collision angle and cost of
%   transport are computed as descriptors, not classifier inputs. Axis convention:
%   1 = medio-lateral, 2 = fore-aft, 3 = vertical.
%
%   See also: GaitSel_DetectGaitEvents, GaitSelMulti_ExportTidyCSV.


    m = B.mass; g = B.g; BW = m*g;
    t = B.time; F = B.force; v = B.comVel; pos = B.com;
    si = E.stepIdx; nStep = numel(si) - 1;
    % Travel direction: the birds cross the runway either way, so fore-aft sign is
    % not consistent across trials. Normalise so forward (direction of travel) is
    % positive, making speed, step length and the acceleration sign comparable.
    travelDir = sign(pos(end,2) - pos(1,2)); if travelDir == 0, travelDir = 1; end

    % --------------------------- per step --------------------------------
    stepRows = {};
    for k = 1:nStep
        ii = si(k):si(k+1);
        if numel(ii) < 4, continue; end
        st = stanceWin(E, k, ii);                    % grounded window for forces
        r = baseId(B);
        r.stepIndex   = k;
        r.stepPeriod  = E.stepPeriod(k);
        r.stepFreq    = 1/E.stepPeriod(k);
        r.groundedFrac= E.groundedFrac(k);
        % Primary per-limb duty factor from foot-marker contact periods
        % (GaitSel_FootContacts): true stance duration / stride period, valid for
        % walking and running alike. Missing or invalid marker contact (a
        % non-finite value, or a zero contact fraction, which is non-physical for a
        % real step and signals a foot-marker detection failure) is recorded as NaN,
        % NOT zero and NOT silently substituted, so those steps drop out of plots
        % and summaries rather than forming a false zero-duty cluster. The
        % force-based grounded-fraction/2 estimate is kept separately in
        % dutyFactor_gf for cross-check.
        r.dutyFactor    = E.dutyMarker(k);
        if ~isfinite(r.dutyFactor) || r.dutyFactor <= 0, r.dutyFactor = NaN; end
        r.dutyFactor_gf = E.groundedFrac(k)/2;
        r.contactTime   = E.contactTimeMarker(k);      % marker-based stance (s)
        r.flightTime  = E.flightTime(k);
        r.hasFlight   = double(E.hasFlight(k));
        r.dutyLimb    = E.dutyLimb(k);
        r.meanSpeed   = mean(v(ii,2))*travelDir;            % forward speed (+)
        r.speed_TD    = v(ii(1),2)*travelDir;
        r.stepLength  = (pos(ii(end),2) - pos(ii(1),2))*travelDir;
        if ~isempty(st)
            r.peakVertForce_BW = max(F(st,3))/BW;
            r.peakForeAftForce_BW = max(abs(F(st,2)))/BW;
            r.peakResultantForce_BW = max(sqrt(sum(F(st,:).^2,2)))/BW;
        end
        r.impulseForeAft = trapz(t(ii), F(ii,2))*travelDir;  % over the step, forward+
        r.accelForeAft   = r.impulseForeAft/(m*E.stepPeriod(k));   % + = speeding up
        r.accelSign      = sign(r.accelForeAft);
        r.fracDV         = r.accelForeAft*E.stepPeriod(k)/r.meanSpeed;   % fractional speed change (fore-aft comparison)
        r.accClass_fa    = classifyAcc(r.fracDV);                        % fore-aft comparison criterion
        % per-step CoM work vs energy change (work-energy identity, short window)
        vv = v(ii,:); z = pos(ii,3); tt = t(ii);
        KE = 0.5*m*sum(vv.^2,2); PE = m*g*z; Pw = sum(F(ii,:).*vv,2);
        r.W_CoM_net = trapz(tt, Pw);
        r.dE_CoM    = (KE(end)+PE(end)) - (KE(1)+PE(1));
        r.dKE_step  = KE(end)-KE(1);  r.dPE_step = PE(end)-PE(1);
        r.fracDE    = r.dE_CoM / max(m*r.meanSpeed^2, eps);  % energy-fraction comparison (previous primary)
        r.fracEG    = r.dE_CoM / max(m*g*r.stepLength, eps); % energy grade dE/(m g L): PRIMARY steadiness metric
        r.accClass  = classifyAcc(r.fracEG, 0.05);           % PRIMARY steadiness (energy grade, |grade| <= 0.05)
        % per-step work-energy identity residual (QC gate): |W_CoM - dE_CoM|
        % relative to the CoM mechanical-energy fluctuation over the step.
        r.WE_relresid = abs(r.W_CoM_net - r.dE_CoM) / max(range(KE+PE), eps);
        % per-step CoM energy exchange (the CoM completes one KE/PE oscillation
        % per step; gaits are symmetric, so the step is the natural unit).
        Ekf = 0.5*m*(vv(:,1).^2 + vv(:,2).^2); Ev = m*g*z + 0.5*m*vv(:,3).^2; Et = Ekf+Ev;
        Wf = posSum(Ekf); Wv = posSum(Ev); Wt = posSum(Et);
        r.recovery  = 100*(Wf+Wv-Wt)/max(Wf+Wv,eps);
        r.congruity = 100*mean(sign(diff(Ekf))==sign(diff(Ev)));
        r.KE_range  = range(Ekf); r.PE_range = range(m*g*z);
        r.ampRatio_KEtoPE = range(Ekf)/max(range(Ev),eps);
        r.vertExcursion   = range(z);
        r.W_CoM_pos = trapz(tt,max(Pw,0)); r.W_CoM_neg = trapz(tt,min(Pw,0));
        % hodograph rotation sense (classifier feature) + collision descriptors
        r.hodoArea = hodographArea(vv, travelDir);
        [r.collisionAngle, r.collisionFraction, r.CoTmech] = ...
            collisionMetrics(F(ii,:), vv, r.meanSpeed, BW);
        [r.gaitObjective, r.mechClass] = classify(r.hasFlight, r.hodoArea);
        [r.legLen_TD, r.legAngle_TD, r.dLegLen, r.dLegAngle] = legGeom(B, ii, st, travelDir);
        stepRows{end+1} = r; %#ok<AGROW>
    end

    % --------------------------- per stride ------------------------------
    starts = 1:2:(nStep-2);                           % non-overlapping tiling: stride s0 = steps s0, s0+1
    strideRows = {};
    for s0 = starts
        ii = si(s0):si(s0+2);                         % stride = two steps
        if numel(ii) < 8, continue; end
        vv = v(ii,:); z = pos(ii,3); tt = t(ii);
        Ekf = 0.5*m*(vv(:,1).^2 + vv(:,2).^2);        % forward kinetic energy
        Ev  = m*g*z + 0.5*m*vv(:,3).^2;               % vertical energy (PE + vert KE)
        Et  = Ekf + Ev;
        Wf = posSum(Ekf); Wv = posSum(Ev); Wt = posSum(Et);
        rec = 100*(Wf+Wv-Wt)/max(Wf+Wv,eps);
        congr = 100*mean(sign(diff(Ekf))==sign(diff(Ev)));
        Pw = sum(F(ii,:).*vv,2);
        sf = s0:s0+1;                                 % the stride's two steps
        r = baseId(B);
        r.strideIndex  = s0;
        r.stridePeriod = tt(end)-tt(1);
        r.strideFreq   = 1/(tt(end)-tt(1));
        r.meanSpeed    = mean(vv(:,2))*travelDir;           % forward speed (+)
        r.strideLength = (pos(ii(end),2)-pos(ii(1),2))*travelDir;  % fore-aft displacement over stride
        r.recovery     = rec;
        r.congruity    = congr;
        r.Wf_horiz = Wf; r.Wv_vert = Wv; r.Wtot_com = Wt;
        r.dKE = Ekf(end)-Ekf(1);  r.dPE = m*g*(z(end)-z(1));
        r.KE_range = range(Ekf);  r.PE_range = range(m*g*z);
        r.ampRatio_KEtoPE = range(Ekf)/max(range(Ev),eps);
        r.vertExcursion = range(z);
        r.W_CoM_pos = trapz(tt,max(Pw,0));  r.W_CoM_neg = trapz(tt,min(Pw,0));
        r.W_CoM_net = trapz(tt, Pw);                        % net external CoM work
        r.dE_CoM    = (Et(end)) - (Et(1));                  % total mechanical energy change
        r.WE_relresid = abs(r.W_CoM_net - r.dE_CoM) / max(Wf+Wv, eps);
        % Steadiness metrics (see GaitSelMulti_QCGate.m / 02_clean.R). PRIMARY = energy grade
        % fracEG = dE_CoM/(m g L): net CoM mechanical energy change per unit distance as a
        % fraction of body weight (speed- and size-neutral effective fore-aft grade). fracDE =
        % dE_CoM/(m v^2) is the previous energy-fraction comparison; fracDV (below) is the
        % fore-aft speed-change comparison.
        r.fracDE = r.dE_CoM / max(m*r.meanSpeed^2, eps);
        r.fracEG = r.dE_CoM / max(m*g*r.strideLength, eps);
        r.accClass = classifyAcc(r.fracEG, 0.05);           % PRIMARY steadiness (energy grade)
        r.hodoArea = hodographArea(vv, travelDir);
        [r.collisionAngle, r.collisionFraction, r.CoTmech] = ...
            collisionMetrics(F(ii,:), vv, r.meanSpeed, BW);
        r.groundedFrac = mean(E.groundedFrac(sf));
        % marker-based per-limb duty factor for the stride (mean of its two steps),
        % falling back to grounded fraction / 2 where markers are missing.
        dm = mean(E.dutyMarker(sf),'omitnan');
        if ~isfinite(dm) || dm <= 0, dm = NaN; end   % invalid marker duty -> NaN, not substituted
        r.dutyFactor    = dm;
        r.dutyFactor_gf = mean(E.groundedFrac(sf))/2;
        r.contactTime   = mean(E.contactTimeMarker(sf),'omitnan');
        r.flightTime   = mean(E.flightTime(sf));
        r.dutyLimb     = mean(E.dutyLimb(sf),'omitnan');
        r.hasFlight    = double(any(E.hasFlight(sf)));
        r.accelForeAft = trapz(tt,F(ii,2))*travelDir/(m*(tt(end)-tt(1)));  % + = speeding up
        r.accelSign    = sign(r.accelForeAft);
        r.fracDV       = r.accelForeAft*r.stridePeriod/r.meanSpeed;        % fractional speed change (fore-aft comparison)
        r.accClass_fa  = classifyAcc(r.fracDV);                            % fore-aft comparison criterion
        [r.gaitObjective, r.mechClass] = classify(r.hasFlight, r.hodoArea);
        r.gaitCodeHand = B.gaitCode;
        strideRows{end+1} = r; %#ok<AGROW>
    end
end

% -------------------------------------------------------------------- helpers
function r = baseId(B)
    r = struct('boutID',B.boutID,'bird',B.bird,'gaitLabelHand',B.gaitLabel, ...
               'trialNo',B.trialNo,'dateCode',B.dateCode,'mass_kg',B.mass, ...
               'driftRMS_vert_mm',1000*B.driftRMS(3),'WE_r',B.WEcheck.r);
end

function st = stanceWin(E, k, ii)
    if isfinite(E.stanceOn(k)) && isfinite(E.stanceOff(k)) && E.stanceOff(k) > E.stanceOn(k)
        st = E.stanceOn(k):E.stanceOff(k);
    else
        st = ii;                                      % fall back to whole step
    end
end

function s = posSum(E), s = sum(max(diff(E),0)); end

function [L0,ang0,dL,dAng] = legGeom(B, ii, st, travelDir)
% Virtual-leg geometry at the start and end of stance: length, angle, and their changes.
%
% The leg ANGLE is measured anticlockwise from the horizontal with fore-aft POSITIVE IN THE
% DIRECTION OF TRAVEL, so it is comparable across trials regardless of which way the bird
% crossed the runway (matching Blum et al. 2014, who report 122.6 deg at touchdown for
% guinea fowl level running). Without the travelDir factor the fore-aft component keeps the
% raw lab sign, which reflects the angle about the vertical for the trials that ran the other
% way and folds the distribution about 90 deg. Leg LENGTH is a magnitude and is unaffected.
    L0=NaN; ang0=NaN; dL=NaN; dAng=NaN;
    if isempty(st), return; end
    if nargin < 4 || ~isfinite(travelDir) || travelDir == 0, travelDir = 1; end
    % contacting foot = the one lower (in vertical) during stance
    zR = mean(B.footR(st,3),'omitnan'); zL = mean(B.footL(st,3),'omitnan');
    if ~isfinite(zR) && ~isfinite(zL); return; end
    if ~isfinite(zL) || (isfinite(zR) && zR <= zL); foot = B.footR; else; foot = B.footL; end
    a = st(1); b = st(end);
    legA = B.com(a,:)-foot(a,:); legB = B.com(b,:)-foot(b,:);
    if all(isfinite(legA)); L0=sqrt(sum(legA.^2)); ang0=atan2d(legA(3),legA(2)*travelDir); end
    if all(isfinite(legB)) && all(isfinite(legA))
        dL=sqrt(sum(legB.^2))-L0; dAng=atan2d(legB(3),legB(2)*travelDir)-ang0; end
end

function [gait,mech] = classify(hasFlight, hodoArea)
% Two threshold-free mechanistic features, no presupposed number of gaits:
%   (1) aerial phase present (hasFlight);
%   (2) hodograph rotation sense = sign of the CoM velocity-loop signed area,
%       positive = pendular/vaulting (walking sense), negative = bouncing
%       (running sense). The boundary is zero (rotation reversal), not a chosen
%       level. The CoM velocity hodograph and the collision view of gait are due
%       to Kuo and the Ruina/Cornell collision-model group (Kuo 2001, 2002, 2007;
%       Ruina, Bertram & Srinivasan 2005; Kuo, Donelan & Ruina 2005; Adamczyk &
%       Kuo 2009); this uses the rotation SENSE as an objective classifier.
% Gaits are the populated combinations: grounded+pendular = walk,
% grounded+bouncing = grounded run, aerial+bouncing = aerial run.
    if hodoArea < 0; mech = "bouncing"; else; mech = "pendular"; end
    if hasFlight;              gait = "aerialRun";
    elseif mech == "bouncing"; gait = "groundedRun";
    else;                      gait = "walk";
    end
end

function a = hodographArea(vv, travelDir)
% Signed area of the CoM velocity loop (shoelace/Green's theorem) over the window.
% x = fore-aft velocity fluctuation about the cycle mean, made positive in the
% direction of travel; y = vertical velocity. Sign is travel-direction invariant.
% Positive = counterclockwise (pendular), negative = clockwise (bouncing).
    x = vv(:,2)*travelDir; x = x - mean(x);
    y = vv(:,3);
    xn = [x(2:end); x(1)]; yn = [y(2:end); y(1)];
    a = 0.5*sum(x.*yn - xn.*y);
end

function [Phi, collFrac, CoTmech] = collisionMetrics(Fw, Vw, meanSpeed, BW)
% Collision-based descriptors (Lee, Comanescu, Butcher & Bertram 2013) from the
% force and velocity vectors over the window (force-magnitude weighting naturally
% downweights the flight phase). theta = deviation of F from vertical; lambda =
% deviation of V from horizontal (travel); phi = deviation of F,V from
% perpendicular. Phi (force-and-velocity-weighted phi) ~ dimensionless mechanical
% cost of transport (CoTmech). Angles use absolute components, so travel-direction
% invariant. Computed over ground-contact where force is non-negligible.
    Fmag = sqrt(sum(Fw.^2,2)); Vmag = sqrt(sum(Vw.^2,2)); dotFV = sum(Fw.*Vw,2);
    ok = Fmag > eps & Vmag > eps;
    theta  = acos(min(1, abs(Fw(:,3))./max(Fmag,eps)));   % F from vertical
    lambda = acos(min(1, abs(Vw(:,2))./max(Vmag,eps)));   % V from horizontal (travel)
    phi    = asin(min(1, abs(dotFV)./max(Fmag.*Vmag,eps)));% F,V from perpendicular
    w = Fmag.*Vmag; w(~ok) = 0;
    Phi      = sum(w.*phi)/max(sum(w),eps);
    collFrac = sum(w.*(phi./max(theta+lambda,eps)))/max(sum(w),eps);
    CoTmech  = mean(abs(dotFV))/max(meanSpeed*BW,eps);     % dimensionless
end

function c = classifyAcc(x, thr)
% Steadiness class from a signed fractional/grade measure: accelerating (x > thr) /
% decelerating (x < -thr) / steady (|x| <= thr; non-finite -> steady). Default thr = 0.10
% (the fore-aft Birn-Jeffery speed-change criterion); pass 0.05 for the energy-grade primary.
    if nargin < 2, thr = 0.10; end
    if ~isfinite(x); c = "steady"; elseif x > thr; c = "accelerating";
    elseif x < -thr; c = "decelerating"; else; c = "steady"; end
end
