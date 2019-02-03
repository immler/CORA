function names = names_of(obj)
    [n, m] = size(obj);
    names = {};
    for i = 1:n
        for j = 1:m
            names = unique([names, obj(i, j).names_of_var]);
        end
    end
end