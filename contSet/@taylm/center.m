function res = center( obj )
    res = arrayfun(@(o)center1(o), obj);
end

function res = center1( obj )
    obj = obj;
    obj.max_order = 0;
    [res, rem] = compress(obj);
    mergeProperties(res, obj, res);
end

