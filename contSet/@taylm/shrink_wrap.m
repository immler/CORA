%
% Reference: Shrink wrapping for Taylor models revisited, Florian Bünger
%            (Numerical Algorithms (2018) 78:1001-1017)
%            https://doi.org/10.1007/s11075-017-0410-1
function res = shrink_wrap(obj)
    [n, m] = size(obj);
    if m ~= 1
        error ('expected column vector tm')
    end
    names = names_of(obj);
    [~, nnames] = size(names);
    if nnames > n
        error ('too many variables - consider extending your taylm')
    end
    if nnames < n
        names = [names, ...
            arrayfun(@(i) {['shrink_wrap_dummy', num2str(i)]}, 1:(n - nnames))]
    end
    % 1. and 2. ================
    [cst, lin, nln] = split_linear(obj)
    p = obj - cst;
    p = set_remainders(p, repmat(interval(0,0), n, 1))
    
    A = zeros(n, n);
    for i = 1:n
        for j = 1:n
            A(i, j) = coeff_of_names(obj(i), names(j), 1);
        end
    end
    A = A
    
    % 3. ================
    if cond(A) > 1e5
        error (['Condition number too high - consider implementing a' ...
            'fallback, e.g., Lemma 3 in Buengers article'])
    end
    
    % 4. ================
    R = inv(A)
    
    % 5. ================
    r = zeros(n,1)
    for i = 1:n
        r(i) = supremum(abs(cst(i).remainder));
    end
    r = r
    rs = abs(R)*r
    
    % 6. ================
    I = arrayfun(@(n) taylm(interval(-1,1), obj(1).max_order, n), names');
    g = R * p - I    
    dg = jacobian(getSyms(g)) % TODO: this can and should be 
    
%    c = center(obj)
%    obj = obj - c
end

function res = verify_shrink(r, rt, s, R)
    if ~(abs(R) * r <= rt)
        disp('scaling of r does not work')
        res = 0
    end
    
end

function names = names_of(obj)
    [n, m] = size(obj);
    names = {};
    for i = 1:n
        for j = 1:m
            names = unique([names, obj(i, j).names_of_var]);
        end
    end
end