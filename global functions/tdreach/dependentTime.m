function [x0, res] = dependentTime(f, samples, h)
    % find closest point
    [n, N] = size(samples); % N sample points of dimension n
    h_eps = h/4;

    x0 = mean(samples, 2);
    nrml = f(x0); % reduce width of set along this direction!
    C = zeros(n * N, n + 1);
    ha = zeros(n + 1, 1);
    A = zeros(2 * N, n + 1);
    b = zeros(2 * N, 1);
    
    for i = 1:N
        xi = samples(:,i);
        dxi = xi - x0;
        Dfn = Jacobian(xi, 1e-5, f) * nrml;
        C((i-1)*n + 1:i*n,1) = Dfn;
        C((i-1)*n + 1:i*n,2:n+1) = Dfn * dxi' + f(xi)*nrml';
        A(i,:) = [1, dxi'];
        A(N + i,:) = -[1, dxi'];
    end
    d = repmat(-nrml, N, 1);
    b = [repmat(h, N, 1); repmat(-h_eps, N, 1)];
    res = lsqlin(C, d, A, b);
end
