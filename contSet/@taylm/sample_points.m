function res = sample_points (obj, N)
% samples
    names = names_of(obj);
    [~, nnames] = size(names);
    ranges = repmat({-1:2/N:1}, size(names));
    outputs = cell(size(ranges));
    [outputs{:}] = ndgrid(ranges{:});
    [n, m] = size(outputs{1});
    res = zeros(nnames, n * m);
    for i = 1:n*m
        sample = zeros(1, nnames);
        for v = 1:nnames
            pv = outputs{v};
            sample(v) = pv(i);
        end
        res(:,i) = point_eval(obj, names, sample);
    end
end