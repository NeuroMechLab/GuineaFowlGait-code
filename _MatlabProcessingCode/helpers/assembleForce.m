function [F, cal] = assembleForce(Fx, Fy, Fz, faKin, t)
%ASSEMBLEFORCE  Auto-calibrate the horizontal force axes to the kinematic travel direction.
%   [F, cal] = assembleForce(Fx, Fy, Fz, faKin, t) picks the horizontal force axis (summed
%   X or Y) whose cumulative impulse best tracks the kinematic fore-aft velocity, assigns it
%   to fore-aft with the sign that makes impulse follow velocity, and puts the other axis on
%   medio-lateral. Because the two rig generations mount their plates differently relative to
%   the motion-capture frame, the plate-to-anatomy axis map is not fixed; matching cumulative
%   horizontal impulse to travel direction recovers it per trial (a per-(study,date) consensus
%   is taken over trials in the batch). faKin is the kinematic fore-aft position; its time
%   derivative is the fore-aft velocity. Returns F = [medio-lateral, fore-aft, vertical] and a
%   cal struct (faAxis "X"/"Y", the two impulse-velocity correlations cx/cy, the winning faCorr,
%   and faSign).
    Fx = Fx(:); Fy = Fy(:); Fz = Fz(:); t = t(:);
    v  = gradient(faKin(:), t);                 % kinematic fore-aft velocity
    Ix = cumtrapz(t, Fx); Iy = cumtrapz(t, Fy);
    cx = corrSafe(Ix, v); cy = corrSafe(Iy, v);
    if abs(cx) >= abs(cy)
        faF = sign(cx)*Fx; mlF = Fy; axisFA = "X"; faCorr = cx;
    else
        faF = sign(cy)*Fy; mlF = Fx; axisFA = "Y"; faCorr = cy;
    end
    F = [mlF, faF, Fz];
    cal = struct('faAxis',axisFA,'cx',cx,'cy',cy,'faCorr',faCorr, ...
                 'faSign',sign(faCorr));
end
