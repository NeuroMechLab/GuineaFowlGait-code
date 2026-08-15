function r = corrSafe(a, b)
%CORRSAFE  Pearson correlation that returns 0 on degenerate input.
%   r = corrSafe(a, b) is corrcoef(a,b) over the jointly finite samples, but returns 0
%   (rather than NaN or erroring) when there are fewer than 3 valid pairs or either input
%   has zero variance. Used by the force-axis auto-calibration (assembleForce).
    a = a(:); b = b(:); ok = isfinite(a) & isfinite(b);
    if nnz(ok) < 3 || std(a(ok)) == 0 || std(b(ok)) == 0, r = 0; return; end
    c = corrcoef(a(ok), b(ok)); r = c(1,2);
end
