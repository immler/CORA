function res = picard_approx(f, init, h, optns)
%    vars = arrayfun (@(i) {['x' num2str(i)]}, 1:dim);
%    vars = {'a' 'b'};
%    init = taylm(interval(-ones(dim,1), ones(dim,1)), optns.picard_order, vars');
    iterate = init;
    error = Inf;
    i = 0;
    while optns.picard_threshold < error
        i = i + 1;
        if mod (i, 10) == 0
            disp(['more than ', num2str(i),' Picard iterations?!'])
        end
        old = iterate;
        iterate = init + integrate_timesubst(f(old), optns.time_var, h);
        error = max(supremum(abs(interval(old - iterate))));
    end
    res = iterate;
end