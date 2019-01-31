% divide taylm by a variable (with index) v, with rest
% obj = r + v * q
function [q, r] = divmod ( obj, v)
    q = obj;
    r = obj;
    [mons, vars] = size(obj.monomials);
    q.remainder = interval(0, 0);
    for i=1:mons
        % part of quotient
        if obj.monomials(i, v) > 0
            r.monomials(i,:) = zeros(1,vars);
            r.coefficients(i) = 0;
            q.monomials(i, 1) = q.monomials(i, 1) - 1; % degree
            q.monomials(i, v) = q.monomials(i, v) - 1; % power of v
        %part of remainder
        else
            q.monomials(i,:) = zeros(1, vars);
            q.coefficients(i) = 0;
        end
    end
    [r, rrem] = compress(r);
    r.remainder = r.remainder + rrem;
    [q, qrem] = compress(q);
    q.remainder = q.remainder + qrem;
end
