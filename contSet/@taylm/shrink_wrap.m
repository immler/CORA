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
    dg = jacobian(getSyms(g)); % TODO: this can and should be done on the level of tms
    dg = arrayfun(@(d) taylm(d, repmat(interval(-1, 1), 1, n)), dg)
    % 7. ================
    small_factor = 1.00000001;
    q_max = 1.01;
    q_tol = 1e-12;
    iter_max = 3;
    q = 1 + rs;
    improve = true;
    iter = 0;
    while improve && iter < iter_max
        s = zeros(size(rs));
        q_old = q;
        for i = 1:n
            for j = 1:n
                % compute upper bound t of |dgi_dxj[-q, q]|
                q_intervals = arrayfun(@(n, qn) taylm(interval(-qn, qn), obj(i).max_order, n), names, q');
                insert = horner(dg(i, j), names, q_intervals);
                insert.max_order = 0;
                [insert, rem] = compress(insert);
                t = supremum(abs(insert.remainder + rem));
                s(i) = s(i) + t * (q(j) - 1);
            end
            q(i) = 1 + rs(i) + s(i);
            if q(i) > q_max
                error (['Shrink wrapping seems not promising (', num2str(q(i)), ' -> try other fallback strategies'])
            end
        end
        improve = any((q-q_old)./q > q_tol);
        iter = iter + 1;
    end
    % 8. ================
    s = small_factor * s;
    q = ones(n,1) + rs + s;
    
    % 9. ================
    
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