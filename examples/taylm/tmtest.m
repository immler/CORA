function tmtest()

    beta=8.0/3.0;
    rho=28;
    sigma=10;
    %res = picard(lorenz, tm);

%    lorenz = [t, (sigma * (y - x)), (x * (rho - z) - y), (x * y - beta * z)];
%    picard([t0 x0 y0 z0], [t0 x0 y0 z0], lorenz, [t x y z], t, 2);
%    lorenz = picardIter([t, x0, y0, z0], lorenz, [t x y z], t, 6)
    
    flow = @(x) [1 + x(1)^2];
    tm = taylm(interval(-1, 1), 4, {'x'});
    res_flow = picardIter(tm, flow, 't', 0.1, 30, 1.5, 1.01)
    
    lorenz = @(x) [sigma * (x(2)-x(1)), x(1) * (rho - x(3)) - x(2), x(1)*x(2) - beta * x(3)];
    tm = taylm(interval([0 0 28], [0.1 0.2 28]), 3, {'x' 'y' 'z'});
    res_lorenz = picardIter(tm, lorenz, 't', 0.009, 30, 1.1, 1.01)
end

% HORNER INSERT (tm for t) TM:
% split constant (w.r.t., t) part and higher part
% horner = new taylm(const)+tm*(higher/t)

function res = picard(init, f, t, h, x)
    % res = arrayfun(@(i, f) collectAndReduce(i + int(subs(f, x, iterate), t, 0, t), [t, init], maxOrder), init, f);
    res = arrayfun(@(i, f) i + integrate(f, t, h), init, f(x));
end

function res = picardIter(init, f, t, h, N, widening_scale, narrowing_scale)
    iterate = init;
    do = true;
    n = 0;
    % widening
    while do
        n = n + 1;
        old = iterate;
        old = set_remainders(old, widening_scale * remainders(old));
        iterate = picard(init, f, t, h, old);
        do = (not (subset(iterate, old)) & n < N);
    end
    widening_iterations = n
    if n == N
        error("Picard iteration not converging: you can try to set a higher widening_scale or lower step size!")
    end
    iterate = old; % this is a certified enclosure (the variant of the previous loop)
    % narrowing
    n = 0;
    do = true;
    while do
        n = n + 1;
        old = iterate;
        iterate = picard(init, f, t, h, old);
        % iterate2 ensures progress by at least a factor of narrowing_scale
        iterate_wide = set_remainders(iterate, narrowing_scale * remainders(iterate));
        do = (subset(iterate_wide, old) & n <= N);
    end
    narrowing_iterations = n
    res = old; % this is a certified enclosure (an invariant of the previous loop)
end