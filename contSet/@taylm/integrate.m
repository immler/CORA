function res = integrate( obj, d, h )
    res = integrate_poly(obj, d);
    time = taylm(sym(d), interval(0, h), obj.max_order);
    [res, time] = rescale_dim(res, time);
    res = horner(res, d, time);
    res.remainder = obj.remainder * interval(0, h);
end
