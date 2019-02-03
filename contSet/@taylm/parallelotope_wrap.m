%
% Reference: Shrink wrapping for Taylor models revisited, Florian Bünger
%            (Numerical Algorithms (2018) 78:1001-1017)
%            https://doi.org/10.1007/s11075-017-0410-1
%
% Lemma 3
function res = parallelotope_wrap(obj, varargin)
    small_factor = 1.001;
    q_max = 1.01;
    q_tol = 1e-8;
    iter_max = 25;
    if ~isempty(varargin)
        error('varargin not supported by ptwrap')
    end
    
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
            arrayfun(@(i) {['shrink_wrap_dummy', num2str(i)]}, 1:(n - nnames))];
    end
    
    [cst, lin, h] = split_linear(obj);
    p = obj - cst;
    p = set_remainders(p, repmat(interval(0,0), n, 1));
    
    A = zeros(n, n);
    for i = 1:n
        for j = 1:n
            A(i, j) = coeff_of_names(obj(i), names(j), 1);
        end
    end
    
    if cond(A) > 1e5
        error (['Condition number ', num2str(cond(A)), ' too high - consider implementing a' ...
            'nother fallback'])
    end
    
    % Optimiziation
    he = set_remainders(h, remainders(cst));
    he = set(he, 'opt_method', 'int');
    prob = inv(A)*he;
    r = supremum(abs(interval(prob)));

    % Finish    
    q = ones(n,1) + r;
    cst = cst.set_remainders(repmat(interval(0, 0), n, 1));
    I = arrayfun(@(n) variable_like(n, obj(1)), names');
    res = horner(lin, names, (q .* I)') + cst;
end