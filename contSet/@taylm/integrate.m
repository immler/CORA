function res = integrate( obj, d, h )
%INTEGRATE Taylor model obj w.r.t dependent variable d and evaluate
%          definite integral at bounds a, b
    % find the common variables 
    res = obj;
    t = taylm(sym(d), interval(-1, 1), obj.max_order);
    set(t, 'opt_method', obj.opt_method);
    [t, res] = rescale_dim(t, res);
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
    res.remainder = res.remainder * interval(0, h);
    obj = obj;
    res = res;
    res = res;
end
