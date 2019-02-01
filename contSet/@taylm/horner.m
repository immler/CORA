%insert the taylor models args into the variables vs of the (vector) taylm obj
function res = horner (obj, vs, args)
    res = arrayfun(@(o) horners(o, vs, args), obj);
end
%insert the taylor models args into the variables vs of the scalar taylm obj
function res = horners ( obj, vs, args )
    [~, sv] = size(vs);
    if sv == 0
        res = compresso(obj);
    else
        w = vs(1, 1);
        ws = vs(1, 2:end);
        brgs = args(1, 2:end);
        i = index_of(obj, w);
        if i == 0
            res = horners(obj, ws, brgs);
        else
            [q, r] = divmod(obj, i);
            r = horners(r, ws, brgs);
            q = horners(q, vs, args);
            a = args(1,1);
            p = a * q;
            res = (r + p);
        end
    end
end

function res = compresso(obj)
    [res, resr] = compress(obj);
    res.remainder = res.remainder + resr;
end

%insert the taylor model arg into the variable v of the scalar taylm obj
% TODO: just legacy code for single variable
function res = s_horner1 ( obj, v, arg )
    i = index_of(obj, v);
    if i == 0
        res = obj;
    else
        [q, r] = divmod(obj, i);
        res = r + arg * s_horner1(q, v, arg);
    end
end
