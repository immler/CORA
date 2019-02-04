function vdp_reach()
    % optns for Picard iteration
    optns.picard_threshold = 1e-6;
    optns.picard_order = 10;
    optns.widening_scale = 1.1;
    optns.narrowing_scale = 1.1;
    optns.time_var = {'t'};
    optns.remainder_estimation = [1e-10; 1e-10];
    optns.parallelotope_factor = 1.1;
    optns.shrinking_mod=1;
    optns.shrink_wrap_options = { 'small_factor' 1.01;
        'iter_max' 5;
        'q_max' 1.01;
        'q_tol' 1e-12};

    % iniitalization of simulation and tdreach
    options.timeStep=0.1;
    options.tFinal=0.1;
    
    tm0 = taylm(interval([1.25; 2.25], [1.55; 2.35]),10, {'x'; 'y'}, 'int', 1e-3, 1e-12);
    zono0 = zono_of_taylm(tm0, ['x', 'y']);
    options.projectedDimensions=[1 2];

    % fixed options for simulation
    options.tStart=0;
    options.tensorOrder = 1;
    options.uTrans = 0;
    options.U = zonotope([0,0]);
    options.x0=center(zono0);
    options.R0=zono0;
    
    mu=1;
    vdpt = @(x) [x(2); mu*(1-x(1)^2)*x(2)-x(1)];

    % experimenting with return time...
    x0 = [1.4; 2.3]
    res = approxReturnTimeDerivative(vdpt, x0, 0.1)
return
    % compute reachable set
    [reach, rs] = timeSeries(tm0, vdpt, options.timeStep, options.tFinal, optns);

    % simulation
    vdp_sys = nonlinearSys(2,1,@vanderPolEq, options);
    %--------------------------------------------------------------------------
    runs = 60;
    fractionVertices = 0.5;
    fractionInputVertices = 0.5;
    inputChanges = 6;
    simRes = simulate_random(vdp_sys, options, runs, fractionVertices, fractionInputVertices, inputChanges);
    
    %plot results--------------------------------------------------------------
    figure;
    hold on
    
    for i=1:rs
        %plot flowpipes
        zono = zono_of_taylm(reach{i}{2}, ['t', 'x', 'y']);
        plotFilled(zono,options.projectedDimensions,[.5 .5 .5],'EdgeColor','black');
    end
    
    %plot initial set
    plotFilled(options.R0,options.projectedDimensions,'w','EdgeColor','black');

    %plot discrete sets
    for i=1:rs
        zono = zono_of_taylm(reach{i}{1}, ['t', 'x', 'y']);
        plotFilled(zono,options.projectedDimensions,[.8 .8 .8],'EdgeColor','black');
    end
    
    %plot discrete as grid
    for i=1:rs
        grid = grid_of_taylm(reach{i}{1}, options.projectedDimensions, 4);
        [sg, ~] = size(grid);
        for j = 1:sg
            zono = zono_of_taylm(grid{j}, ['t', 'x', 'y']);
            plotFilled(zono,[1 2], [1.0 0.5 1.0],'EdgeColor','none');
        end
    end
    
    %plot simulation results      
    for i=1:length(simRes.t)
        plot(simRes.x{i}(:,options.projectedDimensions(1)),simRes.x{i}(:,options.projectedDimensions(2)),'Color',0*[1 1 1]);
    end

end
