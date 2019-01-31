function res = zono_of_taylm( obj, names )
    obj = order_1_taylm(obj)
    [so, ~] = size(obj);
    [~, sn] = size(names);
    zono = zeros(so, sn + 1 + so(1));
    for i = 1:size(obj)
        % center
        zono(i, 1) = coeff_of_names(obj(i), [], []);
        % generators
        for n = 1:sn
            zono(i, n + 1) = coeff_of_names(obj(i), [names(1,n)], [1]);
        end
        % interval remainder
        rem = obj(i).remainder;
        supr = supremum(rem);
        infm = infimum(rem);
        zono(i, 1) = zono(i, 1) + (supr + infm) / 2;
        zono(i, sn + i + 1) = (supr - infm) / 2;
    end
    zono = zono;
    res = zonotope(zono);
end

function res = order_1_taylm(obj)
    res = arrayfun(@s_order_1_taylm, obj);
end

function res = s_order_1_taylm( obj )
    objz = obj;
    objz.max_order = 1;
    [res, rem] = compress(objz);
    res.remainder = res.remainder + rem;
end

