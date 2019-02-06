function res = dependentTime(f, samples, h)
    % find closest point
    [n, N] = size(samples); % N sample points of dimension n
    h_eps = h/4;

    x0 = mean(samples, 2);
    nrml = f(x0); % reduce width of set along this direction!
    A = zeros(2 * N, n + 2);
    b = zeros(2 * N, 1);
    C = zeros(N, n + 2);
    d = zeros(N, 1);

    for i = 1:N
        xi = samples(:,i);
        A(i,2:n+2) = [1, xi'];
        b(i) = h;
        A(N + i,2:n+2) = -[1, xi'];
        b(N + i) = -h_eps;
        C(i,1) = -1;
        C(i,2:n+2) = nrml'*f(xi)*[1, xi'];
        d(i) = -nrml'*xi;
    end
    sol = lsqlin(C, d, A, b);
    res = sol(2:end);
end
