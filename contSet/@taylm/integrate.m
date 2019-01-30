function res = integrate( obj, d, h )
%INTEGRATE Taylor model obj w.r.t dependent variable d and evaluate
%          definite integral at bounds a, b
    % find the common variables 
    res = obj;
    t = taylm(sym(d), interval(0, h), obj.max_order);
    [t, res] = rescale_dim(t, res);
    res = res * t; % this also scales the remainder interval correctly
    t_index = find([res.names_of_var{:}] == d) + 1; % because monomials has as first column the degree
    if isempty(t_index)
        res = res; % this means that res was zero, therefore no time variable
    else
        for i = 1:size(res.coefficients)
            deg_t = res.monomials(i, t_index);
            res.coefficients(i) = res.coefficients(i) / (deg_t + 1);
        end
    end
    % Addition
    %res.coefficients = [factor1.coefficients(:); -factor2.coefficients(:)];
    %res.monomials = [factor1.monomials; factor2.monomials];

    % Merge the properties of the two taylor models
    %res = mergeProperties(res,factor1,factor2);

    % Reduce number of terms of the resulting Taylor model
    %[res] = compress(res);

    %res.remainder = factor1.remainder - factor2.remainder;

end

