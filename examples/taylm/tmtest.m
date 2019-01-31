function tmtest()
    optns.picard_iterations = 30
    optns.widening_scale = 5
    optns.narrowing_scale = 1.1
    optns.time_var = 't'
    
    %res = picard(lorenz, tm);
    tm = taylm(sym(@(x, y) 42 + x^2 + 3 * x * y + 4 * y + x*y^3), interval([-1 -1], [1 1]));
    tm.remainder = interval(1, 2)
    res = horner(tm, 'x', taylm(sym(@(y) y), interval(-1, 1)));
    res = horner(tm, ['y', 'x'], [taylm(sym(@(t) t), interval(-1, 1)), taylm(sym(@(s) s), interval(-1, 1))]);
    time = taylm(sym(@(t) t), interval(-1, 1));
%    lorenz = [t, (sigma * (y - x)), (x * (rho - z) - y), (x * y - beta * z)];
%    picard([t0 x0 y0 z0], [t0 x0 y0 z0], lorenz, [t x y z], t, 2);
%    lorenz = picardIter([t, x0, y0, z0], lorenz, [t x y z], t, 6)
    
    flow = @(x) [1 + x(1)^2];
    tm = taylm(interval(-1, 1), 4, {'x'});
    res_flow = picardIter(tm, flow, 0.04, optns);
    
    beta=8.0/3.0;
    rho=28;
    sigma=10;
    lorenzt = @(x) [sigma * (x(2)-x(1)); x(1) * (rho - x(3)) - x(2); x(1)*x(2) - beta * x(3)];
    tm = taylm(interval([-0.1; -0.1; 28], [0.1; 0.1; 28]),4, {'x'; 'y'; 'z'}, 'int');
%    res_lorenz = picardIter(tm, lorenzt, 0.01, optns)
    [reach, rs] = timeSeries(tm, lorenzt, 0.01, 0.15, optns)

    %set options --------------------------------------------------------------
    options.tStart=0; %start time
    options.tFinal=0.15; %final time
    options.x0=[0; 0; 28]; %initial state for simulation
    options.R0=zonotope([options.x0,[0.1 0 0; 0 0.1 0; 0 0 0]]); %initial state for reachability analysis

    options.timeStep=0.1; %time step size for reachable set computation
    options.plotType='frame';
    options.projectedDimensions=[1 3];
    options.tensorOrder = 1;

    options.timeStep=0.02; %time step size for reachable set computation

    %--------------------------------------------------------------------------


    %obtain uncertain inputs
    options.uTrans = 0;
    options.U = zonotope([0,0]); %input for reachability analysis

    %specify continuous dynamics-----------------------------------------------
    lorenz_sys = nonlinearSys(3,1,@lorenz, options); %initialize tank system
    %--------------------------------------------------------------------------

    %create random simulations; RRTs would provide better results, but are
    %computationally more demanding
    runs = 60;
    fractionVertices = 0.5;
    fractionInputVertices = 0.5;
    inputChanges = 6;
    simRes = simulate_random(lorenz_sys, options, runs, fractionVertices, fractionInputVertices, inputChanges);

    
    %plot results--------------------------------------------------------------
    projectedDimensions=[1 3];
    plotOrder = 20;
    figure;
    hold on
    %plot initial set
    plotFilled(options.R0,projectedDimensions,'w','EdgeColor','k');
    
    for i=1:rs
        %plot flowpipes
        zono = zono_of_taylm(reach{i}{2}, ['t', 'x', 'y', 'z']);
        plotFilled(zono,projectedDimensions,[.5 .5 .5],'EdgeColor','none');
    end
    
    %plot discrete sets
    for i=1:rs
        zono = zono_of_taylm(reach{i}{1}, ['t', 'x', 'y', 'z']);
        plotFilled(zono,projectedDimensions,[.8 .8 .8],'EdgeColor','none');
    end
    
    %plot simulation results      
    for i=1:length(simRes.t)
        plot(simRes.x{i}(:,projectedDimensions(1)),simRes.x{i}(:,projectedDimensions(2)),'Color',0*[1 1 1]);
    end

end

function [reach, rs] = timeSeries(x, f, h, T, optns)
    N = T / h + 1;
    reach = cell(N, 2);
    t = 0.0
    i = 1
    while t < T
        dh = min(h, T - t);
        [flowpipe, x] = timeStep(x, f, dh, optns);
        reach{i}{1} = x;
        reach{i}{2} = flowpipe;
        t = t + dh
        i = i + 1;
    end
    rs = i - 1
end

function [flowpipe, final] = timeStep(x, f, h, optns)
    flowpipe = picardIter(x, f, h, optns);
    timestep = taylm(interval(h, h));
    final = horner(flowpipe, 't', timestep);
end

function res = picard(init, f, t, h, x)
    % res = arrayfun(@(i, f) collectAndReduce(i + int(subs(f, x, iterate), t, 0, t), [t, init], maxOrder), init, f);
    x = x;
    fx = f(x);
    res = arrayfun(@(i, f) i + integrate(f, t, h), init, f(x));
end

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