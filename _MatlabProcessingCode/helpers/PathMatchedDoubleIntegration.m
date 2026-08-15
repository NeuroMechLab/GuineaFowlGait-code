function out = PathMatchedDoubleIntegration(force, time, bodyMass, comProxy, opts)
%PATHMATCHEDDOUBLEINTEGRATION  CoM trajectory from force plate, path-matched to kinematics.
%
%   out = PathMatchedDoubleIntegration(force, time, bodyMass, comProxy)
%   out = PathMatchedDoubleIntegration(force, time, bodyMass, comProxy, opts)
%
%   Reconstructs the body center-of-mass (CoM) trajectory from ground reaction
%   force by Newton's second law and double integration. Double integration of
%   noisy force amplifies any baseline error into unbounded position drift, so
%   the integration constants are not assumed: for each axis the initial
%   position, initial velocity, and a constant acceleration-baseline offset are
%   optimized (Nelder-Mead, fminsearch) to MINIMIZE the drift between the
%   force-derived path and an independent kinematics-derived CoM path
%   ("path matching"). This is the lab's standard force-plate CoM reconstruction
%   and a validation primitive (the residual drift measures force/kinematics
%   consistency).
%
%   PHYSICS AND SIGN CONVENTION
%     a(t)   = force(t) / bodyMass + gravityVector        (F = m a)
%     v(t)   = v0 + integral a dt                          (cumulative trapezoid)
%     s(t)   = s0 + integral v dt
%   Ground reaction force is positive UP; gravity is a VECTOR pointing DOWN
%   (-g on the vertical axis, 0 on the others). So a supporting force of m*g at
%   quiet standing gives a = (+m g)/m + (-g) = 0. The acceleration offset absorbs
%   a small residual force-plate baseline error (a(t) becomes a(t) + offset
%   before integrating).
%
%   INPUTS
%     force    : [n x d] ground reaction force (N), one column per axis.
%     time     : [n x 1] time stamps (s); may be non-uniform.
%     bodyMass : scalar body mass (kg).
%     comProxy : [n x d] kinematics-derived CoM position (m) used as the
%                path-match target (e.g. a marker-based CoM estimate). Same size
%                as force. This is what anchors the integration constants.
%     opts     : (optional) struct. Because experiments differ in which lab axis
%                is vertical, you MUST tell the tool how gravity aligns — give
%                either verticalAxis or an explicit gravityVector:
%        .verticalAxis   index (1..d) of the vertical column. The gravity vector
%                        is built as zeros with -g on this axis (pointing down).
%        .g              gravitational acceleration MAGNITUDE (default 9.81).
%                        Used with verticalAxis to build the vector as -g.
%        .gravityVector  [1 x d] explicit gravity vector (override), for oblique
%                        axes or full control. Must point down: its vertical
%                        entry is negative (e.g. [0 0 -9.81]).
%        .optimizeOffset logical, optimize the acceleration baseline offset
%                        (default true). If false the offset is fixed at 0 and
%                        only [s0 v0] are optimized.
%        .initialVel     [1 x d] initial-velocity guess for the optimizer
%                        (default: slope of a line fit to the first few comProxy
%                        samples per axis).
%        .plot           logical, draw the proxy-vs-reconstructed diagnostic
%                        figure (default false). Handle returned in out.fig.
%        .axisLabels     1 x d cellstr for the plot y-labels
%                        (default {'axis 1', 'axis 2', ...}).
%
%   OUTPUT (struct)
%     out.position          [n x d] reconstructed CoM position (m)
%     out.velocity          [n x d] reconstructed CoM velocity (m/s)
%     out.acceleration      [n x d] CoM acceleration used (m/s^2, offset applied)
%     out.initialConditions [3 x d] optimized [s0; v0; accelOffset] per axis
%     out.SSE               [1 x d] sum of squared position drift per axis
%     out.rmseDrift         [1 x d] RMS position drift vs comProxy per axis (m)
%     out.maxDrift          [1 x d] max absolute position drift per axis (m)
%     out.maxDriftIdx       [1 x d] sample index of the max drift per axis
%     out.fig               figure handle if opts.plot is true, else []
%
%   Requires only base MATLAB (fminsearch, cumtrapz, optimset, polyfit).
%
%   Neuromech Lab kit — source: Goldsmith, Hall & Daley 2026 (GF turning),
%   TurnExperiments_BatchCalcTrajectory.m internal DoubleIntegration (renamed
%   here for clarity; struct coupling, modal dialog and forced PDF output
%   removed; generalized to d axes). CC0 1.0.

if nargin < 5, opts = struct(); end
[n, d] = size(force);
if ~isequal(size(comProxy), [n d])
    error('PathMatchedDoubleIntegration:size', 'comProxy must be the same size as force ([n x d]).');
end
if numel(time) ~= n
    error('PathMatchedDoubleIntegration:time', 'time must have n = size(force,1) elements.');
end
if any(~isfinite(force(:))) || any(~isfinite(comProxy(:)))
    error('PathMatchedDoubleIntegration:nan', ...
          'force and comProxy must be finite (trim NaN padding to a clean window first).');
end
time = time(:);

gravityVector  = resolveGravity(opts, d);         % [1 x d], points down on vertical
optimizeOffset = getOpt(opts, 'optimizeOffset', true);
doPlot         = getOpt(opts, 'plot',           false);
axisLabels     = getOpt(opts, 'axisLabels', arrayfun(@(k) sprintf('axis %d', k), 1:d, 'UniformOutput', false));

% F = m a (per axis), including the gravitational acceleration term
accelForce = force ./ bodyMass + gravityVector;   % implicit expansion, [n x d]

% Initial-condition guesses: position and velocity from the kinematic proxy
initVelGuess = getOpt(opts, 'initialVel', startingSlope(comProxy, time));
ICguess = [comProxy(1, :); initVelGuess(:)'; zeros(1, d)];   % [3 x d] = [s0; v0; offset]

% Optimize ICs per axis to minimize drift from the kinematic proxy
finalICs      = nan(3, d);
% Tight convergence: the objective has only 3 free parameters (s0, v0, offset)
% per axis, so it cannot overfit; tightening from fminsearch defaults just finds
% the genuine best-match initial conditions even from a far starting guess.
fminoptions   = optimset('FunValCheck', 'on', 'TolFun', 1e-12, 'TolX', 1e-12, ...
                         'MaxFunEvals', 1e4, 'MaxIter', 1e4);
for k = 1:d
    if optimizeOffset
        obj  = @(x) pmObjective([x(1) x(2) x(3)], accelForce(:,k), comProxy(:,k), time);
        x0   = ICguess(:,k)';
        newX = fminsearch(obj, x0, fminoptions);
        finalICs(:,k) = newX(:);
    else
        obj  = @(x) pmObjective([x(1) x(2) 0], accelForce(:,k), comProxy(:,k), time);
        x0   = ICguess(1:2,k)';
        newX = fminsearch(obj, x0, fminoptions);
        finalICs(:,k) = [newX(:); 0];
    end
end

% Reconstruct the CoM trajectory with the optimized initial conditions
position     = nan(n, d);
velocity     = nan(n, d);
acceleration = nan(n, d);
SSE          = nan(1, d);
rmseDrift    = nan(1, d);
maxDrift     = nan(1, d);
maxDriftIdx  = nan(1, d);

for k = 1:d
    acceleration(:,k)          = accelForce(:,k) + finalICs(3,k);
    [position(:,k), velocity(:,k)] = integrateTwice(acceleration(:,k), finalICs(1:2,k), time);

    drift          = position(:,k) - comProxy(:,k);
    SSE(k)         = sum(drift.^2);
    rmseDrift(k)   = sqrt(mean(drift.^2));
    [maxDrift(k), maxDriftIdx(k)] = max(abs(drift));
end

out = struct('position', position, 'velocity', velocity, 'acceleration', acceleration, ...
             'initialConditions', finalICs, 'SSE', SSE, 'rmseDrift', rmseDrift, ...
             'maxDrift', maxDrift, 'maxDriftIdx', maxDriftIdx, 'fig', []);

if doPlot
    out.fig = plotPathMatch(time, comProxy, position, maxDriftIdx, maxDrift, axisLabels);
end
end


% ======================================================================= local
function fval = pmObjective(ICs, accelForce, proxy, time)
% fminsearch objective: SSE between the force-integrated path and the kinematic
% proxy. ICs = [initialPosition, initialVelocity, accelerationOffset].
adjAccel  = accelForce + ICs(3);                         % adjust acceleration baseline
posForce  = integrateTwice(adjAccel, ICs(1:2), time);    % integrate to position
fval      = sum((posForce - proxy).^2);                  % sum of squared error
end


function [position, velocity] = integrateTwice(acceleration, ICs, time)
% Cumulative-trapezoid double integration. ICs = [initialPosition, initialVelocity].
velocity = ICs(2) + cumtrapz(time, acceleration);
position = ICs(1) + cumtrapz(time, velocity);
end


function v0 = startingSlope(comProxy, time)
% Robust initial-velocity guess: slope of a straight-line fit to the first few
% samples of each proxy column (only a starting point for the optimizer).
d  = size(comProxy, 2);
kk = min(5, size(comProxy, 1));
v0 = zeros(1, d);
for k = 1:d
    p = polyfit(time(1:kk), comProxy(1:kk, k), 1);
    v0(k) = p(1);
end
end


function val = getOpt(opts, name, default)
if isstruct(opts) && isfield(opts, name) && ~isempty(opts.(name))
    val = opts.(name);
else
    val = default;
end
end


function gvec = resolveGravity(opts, d)
% Build the downward gravity VECTOR aligned to the caller's vertical axis.
% Priority: explicit gravityVector > (verticalAxis + g magnitude). Never guesses.
if isstruct(opts) && isfield(opts, 'gravityVector') && ~isempty(opts.gravityVector)
    gvec = opts.gravityVector(:)';
    if numel(gvec) ~= d
        error('PathMatchedDoubleIntegration:gravity', ...
              'gravityVector must have d = %d elements (one per axis column).', d);
    end
    return;
end

g = getOpt(opts, 'g', 9.81);
if ~isscalar(g) || ~(g > 0)
    error('PathMatchedDoubleIntegration:g', ...
          'g must be a positive scalar magnitude; the vector is built as -g on the vertical axis.');
end

if isstruct(opts) && isfield(opts, 'verticalAxis') && ~isempty(opts.verticalAxis)
    va = opts.verticalAxis;
    if ~isscalar(va) || va < 1 || va > d || va ~= round(va)
        error('PathMatchedDoubleIntegration:verticalAxis', ...
              'verticalAxis must be an integer column index in 1..%d.', d);
    end
    gvec = zeros(1, d);
    gvec(va) = -g;            % gravity points down
    return;
end

error('PathMatchedDoubleIntegration:axes', ...
      ['Specify which axis is vertical so gravity aligns: set opts.verticalAxis ' ...
       '(the vertical column index) or give opts.gravityVector explicitly. ' ...
       'The tool does not guess the vertical axis.']);
end


function fh = plotPathMatch(time, proxy, position, maxDriftIdx, maxDrift, axisLabels)
% Diagnostic: kinematic proxy vs. reconstructed CoM per axis, with the max-drift
% instant marked. No modal dialog, no file output (caller decides what to save).
d  = size(proxy, 2);
fh = figure('Color', 'w');
for k = 1:d
    subplot(d, 1, k); hold on;
    hP = plot(time, proxy(:,k),    'k', 'LineWidth', 1.5);
    hR = plot(time, position(:,k), 'r', 'LineWidth', 1.5);
    mi = maxDriftIdx(k);
    plot([time(mi) time(mi)], [proxy(mi,k) position(mi,k)], 'b', 'LineWidth', 2);
    ylabel(axisLabels{k});
    if k == 1, title('Path-matched CoM reconstruction'); end
    if k == d, xlabel('Time (s)'); end
    legend([hP hR], 'kinematic proxy', ...
           sprintf('reconstructed (max drift %.0f mm)', maxDrift(k)*1000), ...
           'Location', 'best'); legend boxoff;
    xlim([time(1) time(end)]);
end
end
