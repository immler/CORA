function lorenz_reach()
    % optns for Picard iteration
    optns.picard_iterations = 30
    optns.widening_scale = 5
    optns.narrowing_scale = 1.1
    optns.time_var = 't'
    
    % iniitalization of simulation and tdreach
    options.tFinal=0.15;
    tm0 = taylm(interval([-0.1; -0.1; 28], [0.1; 0.1; 28]),4, {'x'; 'y'; 'z'}, 'int');
    zono0 = zono_of_taylm(tm0, ['x', 'y', 'z']);
    options.projectedDimensions=[1 3];

    % fixed options for simulation
    options.tStart=0;
    options.tensorOrder = 1;
    options.uTrans = 0;
    options.U = zonotope([0,0]);
    options.x0=center(zono0);
    options.R0=zono0;
    
    
    beta=8.0/3.0;
    rho=28;
    sigma=10;
    lorenzt = @(x) [sigma * (x(2)-x(1)); x(1) * (rho - x(3)) - x(2); x(1)*x(2) - beta * x(3)];
    
    % compute reachable set
    [reach, rs] = timeSeries(tm0, lorenzt, 0.01, 0.01, optns);


    % simulation
    lorenz_sys = nonlinearSys(3,1,@lorenz, options); %initialize tank system
    %--------------------------------------------------------------------------
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
