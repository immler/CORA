% given x(s) with -1 <= s <= 1
%
% the idea:
% s (with standard domain [-1, 1]) encodes time in a Picard iteration on an
% interval [0; h]
% therefore we assume the substitution:
% s = 2*t/h - 1
% t = h * (s + 1) / 2
%
% this function integrates x(s) = x(2*t/h - 1) from 0 to t assuming 0 <= t <= h
% (but always keeps representing time with s)
function res = integrate_timesubst( xs, t, h)
    res = arrayfun(@(x) s_integrate( x, t, h), xs);
end

function res = s_integrate( x, t, h )
    integrand = h/2 * x;
    prim = primitive(integrand, t);
    arg = taylm(interval(-1,-1), x.max_order);
    lower = horner(prim, t, arg);
    res = prim - lower;
    res.remainder = integrand.remainder * 2;
end
