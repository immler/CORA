function res = order_1(obj)
    res = arrayfun(@s_order_1_taylm, obj);
end

function res = s_order_1_taylm( obj )
    objz = obj;
    objz.max_order = 1;
    [res, rem] = compress(objz);
    res.remainder = res.remainder + rem;
    res = merge_in(res, obj);
end

