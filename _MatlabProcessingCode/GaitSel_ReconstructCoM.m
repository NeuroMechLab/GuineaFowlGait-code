function B = GaitSel_ReconstructCoM(B, opts)
%GAITSEL_RECONSTRUCTCOM  Whole-bout path-matched CoM reconstruction.
%
%   B = GaitSel_ReconstructCoM(B) reconstructs the center-of-mass trajectory for
%   an imported bout struct B (from GaitSelMulti_ImportBout) by F = m a and double
%   integration, with the integration constants path-matched to the marker CoM
%   proxy over the WHOLE bout at once (one [pos0; vel0; accelOffset] per axis,
%   not per stride). Vertical is axis 3.
%
%   opts (optional):
%     .g            gravity magnitude (default 9.81).
%
%   Adds to B: com, comVel, comAcc [n x 3]; comICs [3 x 3]; driftRMS, driftMax
%   [1 x 3] (force/kinematics consistency); and the work-energy check WEcheck.
%
%   See also: PathMatchedDoubleIntegration, GaitSel_PerStepStrideMeasures.

    if nargin < 2, opts = struct(); end
    g = getOpt(opts,'g',9.81);

    rec = PathMatchedDoubleIntegration(B.force, B.time, B.mass, B.comProxy, ...
              struct('verticalAxis',3,'g',g));

    B.com     = rec.position;
    B.comVel  = rec.velocity;
    B.comAcc  = rec.acceleration;
    B.comICs  = rec.initialConditions;
    B.driftRMS= rec.rmseDrift;
    B.driftMax= rec.maxDrift;
    B.g = g;

    % work-energy consistency over the whole bout: d(KE+PE) vs integral(F.v)
    KE = 0.5*B.mass*sum(B.comVel.^2, 2);
    PE = B.mass*g*B.com(:,3);
    Pw = sum(B.force .* B.comVel, 2);                 % CoM power from GRF (W)
    Wcum = cumtrapz(B.time, Pw);
    dE   = (KE+PE) - (KE(1)+PE(1));
    ok = isfinite(Wcum) & isfinite(dE);
    r = corr(Wcum(ok), dE(ok));
    B.WEcheck = struct('r', r, 'W_end', Wcum(end), 'dE_end', dE(end));
end

function v = getOpt(o,n,d)
    if isstruct(o) && isfield(o,n) && ~isempty(o.(n)); v = o.(n); else; v = d; end
end
