function res = index_of (obj, n)
    res = find(ismember(obj.names_of_var, n));
    if isempty(res)
        res = 0;
    else
         res = res +1; % +1 because monomials has as first column the degree
    end
end
