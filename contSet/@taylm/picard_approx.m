function res = picard_approx(f, init, h, optns)
%    vars = arrayfun (@(i) {['x' num2str(i)]}, 1:dim);
%    vars = {'a' 'b'};
%    init = taylm(interval(-ones(dim,1), ones(dim,1)), optns.picard_order, vars');
    % do the Picard iteration without remainders...
    init = set_remainders(init, arrayfun(@(o) interval(0,0), init));    
    iterate = init;
    error = Inf;
    i = 0;
    while optns.picard_threshold < error
        i = i + 1;
        if mod (i, 10) == 0
            disp(['more than ', num2str(i),' Picard iterations?!'])
            iterate;
            oldz;
            asdf = oldz-iterate;
            error;
        end
        old = iterate;
        iterate = init + integrate_timesubst(f(old), optns.time_var, h);
        oldz = set_remainders(old, arrayfun(@(o) interval(0,0), old));
        iteratez = set_remainders(iterate, arrayfun(@(o) interval(0,0), iterate));
        error = max(supremum(abs(interval(oldz - iteratez))));
    end
    res = iterate;
end