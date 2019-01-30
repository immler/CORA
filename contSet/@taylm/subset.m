function res = subset(obj1, obj2)
    res = all(arrayfun(@(o1, o2) subset1(o1, o2), obj1, obj2));
end

function res = subset1(obj1, obj2)
    if (size(obj1.coefficients) == size(obj2.coefficients) & ...
            size(obj1.monomials) == size(obj2.monomials))
        if (obj1.coefficients == obj2.coefficients & ...
           obj1.monomials == obj2.monomials & ...
           obj1.remainder <= obj2.remainder)
            res = 1.0;
        else
            res = 0.0;
        end
    else
        res = 0.0;
    end
end
