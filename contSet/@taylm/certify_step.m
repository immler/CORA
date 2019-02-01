function res = certify_step(f, init, h, optns)
    picard = picard_approx(f, init, h, optns);

    % remainder estimation:
    rem = arrayfun(@(e) interval(-e, e), optns.remainder_estimation);
    %rem = [interval(-5.09307e-5, 7.86167e-5); interval(-1.75707e-4, 1.60933e-4)];

    % widening phase:
    i = 0;
    fixed_point = 0;
    while ~fixed_point
        i = i + 1;
        if mod (i, 10) == 0
            disp(['more than ', num2str(i),' widening iterations?!'])
        end
        cand = set_remainders(picard, rem);
        picard2 = init + integrate_timesubst(f(cand), optns.time_var, h);
        
        % subtract plain (!) polynomial, without remainder!
        % really?
        polydiff = picard2 - picard;
        diffi = interval(polydiff);
        if diffi <= rem
            fixed_point = 1;
        else
            rem = optns.widening_scale * diffi;
        end
    end
    
    % narrowing phase:
    i = 0;
    do = 1;
    while do
        i = i + 1;
        if mod (i, 10) == 0
            disp(['more than ', num2str(i),' narrowing iterations?!'])
        end
        cand = set_remainders(picard, rem);
        picard2 = init + integrate_timesubst(f(cand), optns.time_var, h);
        
        % subtract plain (!) polynomial, without remainder!
        polydiff = interval(picard2 - picard);
        diffi = max (1, optns.narrowing_scale) * polydiff;
        if diffi <= rem
            rem = polydiff;
        else
            do = false;
        end
    end
    res = cand;    
end
