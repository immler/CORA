function res = set_remainders(obj, remainders)
%SET_REMAINDERS Summary of this function goes here
%   Detailed explanation goes here
    res = obj;
    for i = 1:size(res)
        res(i).remainder = remainders(i);
    end
end

