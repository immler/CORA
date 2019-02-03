function v = variable_like(n, obj)
    v = merge_in(taylm(interval(-1, 1), 0, n), obj);
end
