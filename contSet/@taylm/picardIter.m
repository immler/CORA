function res = picardIter(init, f, h, optns)
    t = optns.time_var;
    N = optns.picard_iterations;
    widening_scale = optns.widening_scale;
    narrowing_scale = optns.narrowing_scale;
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

function res = picard(init, f, t, h, x)
    % res = arrayfun(@(i, f) collectAndReduce(i + int(subs(f, x, iterate), t, 0, t), [t, init], maxOrder), init, f);
    x = x;
    fx = f(x);
    res = arrayfun(@(i, f) i + integrate(f, t, h), init, f(x));
end
