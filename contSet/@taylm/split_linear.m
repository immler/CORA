function [ const, linear, nonlinear] = split_linear(obj)
    [n, w] = size(obj);
    if w ~= 1
        error ('split_linear on non-column-vector taylor model')
    end
    ivl = taylm(interval(0,0));
    const = repmat(ivl, n, 1);
    linear = repmat(ivl, n, 1);
    nonlinear = repmat(ivl, n, 1);
    for i = 1:n
        [const1, linear1, nonlinear1] = s_split_linear(obj(i, 1));
        const(i, 1) = const1;
        linear(i, 1) = linear1;
        nonlinear(i, 1) = nonlinear1;
    end
end

function [c, l, n] = s_split_linear ( obj )
    c = obj;
    l = obj;
    n = obj;
    [mons, vars] = size(obj.monomials);
    l.remainder = interval(0, 0);
    n.remainder = interval(0, 0);
    for i=1:mons
        if ~(obj.monomials(i, 1) == 0)
            c.monomials(i,:) = zeros(1,vars);
            c.coefficients(i) = 0;
        end
        if ~(obj.monomials(i, 1) == 1)
            l.monomials(i,:) = zeros(1,vars);
            l.coefficients(i) = 0;
        end
        if obj.monomials(i, 1) <= 1
            n.monomials(i,:) = zeros(1,vars);
            n.coefficients(i) = 0;
        end
    end
    [c, rem] = compress(c);
    c.remainder = c.remainder + rem;
    [l, rem] = compress(l);
    l.remainder = l.remainder + rem;
    [n, rem] = compress(n);
    n.remainder = n.remainder + rem;
end
