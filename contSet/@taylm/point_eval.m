function res = point_eval( obj, vars, values )
	res = arrayfun(@(a) s_point_eval(a, vars, values), obj, 'UniformOutput', 0);
    A = cat(1, res{:});
    res = reshape(A, size(res));
end

function res = s_point_eval( obj, vars, values )
    
    % get coefficients
    c = obj.coefficients;
    % get monomials
    degs = obj.monomials(:, 2:end);
    % get var names
    names = obj.names_of_var;
    % get valuation of variables
    v = arrayfun(@(n) values(indexmember(n, vars)), names);
    % make a syms expression
    res = sum(c.*prod(repmat(v,[size(c,1) 1]).^degs,2));
    
end

function res = indexmember(x, xs)
    [~, res] = ismember(x, xs);
end