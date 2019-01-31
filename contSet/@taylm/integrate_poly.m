% iteration of plain polynomial
function res = integrate_poly( obj, d )
    t = taylm(sym(d), interval(-1, 1), obj.max_order);
    [time, res] = rescale_dim(t, obj);
    t_index = find([res.names_of_var{:}] == d) +1; % +1 because monomials has as first column the degree
%    t_index = find(cellfun(@(x)isequal(x,d),res.names_of_var)) + 1; % +1 because monomials has as first column the degree
    if isempty(t_index)
        res = res; % this means that res was zero, therefore no time variable
    else
        for i = 1:size(res.coefficients)
            deg_t = res.monomials(i, t_index);
            res.monomials(i, t_index) = deg_t + 1;
            res.coefficients(i) = res.coefficients(i) / (deg_t + 1);
        end
    end
    res.remainder=interval(0,0);
end
