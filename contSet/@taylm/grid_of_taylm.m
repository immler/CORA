function res = grid_of_taylm(obj, projdims, max_splits )
    p1 = projdims(1);
    p2 = projdims(2);
    names = unique([obj(p1).names_of_var, obj(p2).names_of_var]);
    [~, n_vars] = size(names);
    n_entries = max_splits^n_vars;
    grid = cell(n_entries, 1);
    w = 2 / max_splits;
    for i = 1:n_entries
        d = (i - 1);
        bnds = cell(1, n_vars);
        for v = 1:n_vars
            k = mod(d, max_splits);
            d = floor(d / (max_splits));
            b = -1.0 + k * w;
            bnds{1, v} = taylm(interval(b, b + w), obj(p1).max_order, names(v));
        end
        arg = cellfun(@(b) b, bnds);
        res = horner(obj, names, arg);
        grid{i} = res;
    end
    res = grid;
end