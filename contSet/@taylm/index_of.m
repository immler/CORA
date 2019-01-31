function res = index_of (obj, n)
    res = find(cellfun(@(x)isequal(x,n),obj.names_of_var));
    if isempty(res)
        res = 0;
    else
         res = res +1; % +1 because monomials has as first column the degree
    end
end
