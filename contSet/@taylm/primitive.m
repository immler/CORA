% primitive (for integration w.r.t. d) of polynomial
function res = primitive( obj, d )
    res = arrayfun(@(o) s_primitive(o, d), obj);
end

function res = s_primitive( obj, d )
    t = taylm(interval(-1, 1), obj.max_order, d);
    [res, ~] = rescale_dim(obj, t);
%    t_index = find([res.names_of_var{:}] == d) +1; % +1 because monomials has as first column the degree
    [mons, ~] = size(res.monomials);
    t_index = index_of(res, d); % +1 because monomials has as first column the degree
    if t_index > 0
        for i = 1:mons
            rci = res.coefficients(i);
            if rci ~= 0
                deg_t = res.monomials(i, t_index);
                res.monomials(i, t_index) = deg_t + 1;
                res.monomials(i, 1) = res.monomials(i, 1) + 1;
                res.coefficients(i) = rci / (deg_t + 1);
            end
        end
    end
    res.remainder=interval(0,0);
end
