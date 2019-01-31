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
    [reach, rs] = timeSeries(tm, lorenzt, 0.01, 0.01, optns)

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
