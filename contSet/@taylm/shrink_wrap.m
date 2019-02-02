%
% Reference: Shrink wrapping for Taylor models revisited, Florian Bünger
%            (Numerical Algorithms (2018) 78:1001-1017)
%            https://doi.org/10.1007/s11075-017-0410-1
function res = shrink_wrap(obj)
    [n, m] = size(obj);
    if m ~= 1
        error ('expected column vector tm')
    end
    maxorder = obj(1).max_order;
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
    I = arrayfun(@(n) variable_like(n, obj(1)), names');
    g = R * p - I    
    dg = jacobian(getSyms(g)); % TODO: this can and should be done on the level of tms
    dg = arrayfun(@(d) merge_in(taylm(d, repmat(interval(-1, 1), 1, n)), obj(1)), dg)
    % 7. ================
    small_factor = 1.001;
    q_max = 1.01;
    q_tol = 1e-8;
    iter_max = 20;
    q = 1 + rs;
    improve = true;
    iter = 0;
    while improve && iter < iter_max
        s = zeros(size(rs));
        q_old = q;
        for i = 1:n
            for j = 1:n
                % compute upper bound t of |dgi_dxj[-q, q]|
                q_intervals = arrayfun(@(n, qn) merge_in(taylm(interval(-qn, qn), 0, n), obj(1)), names, q');
                insert = horner(dg(i, j), names, q_intervals);
                insert.max_order = 0;
                insert = insert
                [insert, rem] = compress(insert);
                t = supremum(abs(insert.remainder + rem))
                s(i) = s(i) + t * (q(j) - 1);
            end
            rsi = rs(i)
            si = s(i)
            q(i) = 1 + rs(i) + s(i);
            if q(i) > q_max
                error (['Shrink wrapping seems not promising (', num2str(q(i)), ' -> try other fallback strategies'])
            end
        end
        improve = any((q-q_old)./q > q_tol);
        iter = iter + 1
    end
    % 8. ================
    s = small_factor * s;
    
    % 9. ================
    q = ones(n,1) + rs + s
    if ~(abs(R) * r <= rs)
        error('scaling of r does not work')
    end
    qB = arrayfun(@(q, n) merge_in(taylm(interval(-q, q), 0, n), obj(1)), q', names);
    adgqB = repmat(interval(0,0), n, n);
    for i = 1:n
        for j = 1:n
            dgqBij = horner(dg(i, j), names, qB);
            dgqBij.max_order = 0;
            [c, rem] = compress(dgqBij);
            adgqB(i, j) = abs(c.remainder + rem);
        end
    end
    b = supremum(adgqB * (q - ones(n, 1)));
    s = s;
    if b <= s
        cst = cst.set_remainders(repmat(interval(0, 0), n, 1));
        res = horner(p, names, (q .* I)') + cst;
    else
        b
        s
        error(' last check did not work out as planned...')
    end
%    c = center(obj)
%    obj = obj - c
end

function res = merge_in(tm, obj)
    res = mergeProperties(tm, obj, tm);
end

function v = variable_like(n, obj)
    v = merge_in(taylm(interval(-1, 1), 0, n), obj);
end

function q = verify_shrink(n, rs, s, R)
    
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