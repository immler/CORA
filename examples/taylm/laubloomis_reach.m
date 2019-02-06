function laubloomis_reach()
    % optns for Picard iteration
    optns.picard_threshold = 1e-8;
    optns.picard_order = 4;
    optns.widening_scale = 1.1;
    optns.narrowing_scale = 1.1;
    optns.time_var = {'t'};
    optns.remainder_estimation = 5e-6 * ones(7,1);
    optns.parallelotope_factor = 1.1;
    optns.shrinking_mod = 10;
    optns.shrink_wrap_options = { 'small_factor' 1.01;
        'iter_max' 5;
        'q_max' 1.01;
        'q_tol' 1e-12};
    % iniitalization of simulation and tdreach
    options.timeStep=0.1;
    options.tFinal=1.5;
    
    x0 = [1.2; 1.05; 1.5; 2.4; 1; 0.1; 0.45];
    W = 0.1;
    ivl0 = x0 + W*interval(-ones(7, 1), ones(7, 1));
    
    [ts, xs, ws] = simulateTimeDepReach(ivl0, 1, @ll, 0.05, 20, 1);
    for i = 1:7
        for j=i+1:6
            figure('Name', ['Figure ', num2str(i), ', ', num2str(j)]);
%            plot(xs{i}(end-127:end), xs{j}(end-127:end), 'x')
            plot(xs{i}, xs{j}, 'x')
        end
    end
    figure('Name', 'widths');
    plot(ts, ws)
    return
    names = {'x1'; 'x2'; 'x3'; 'x4'; 'x5'; 'x6'; 'x7'};
    tm0 = taylm(ivl0, 4, names, 'int', 1e-3, 1e-12);
    zono0 = zono_of_taylm(tm0, names');
    options.projectedDimensions=[1 4];

    % fixed options for simulation
    options.tStart=0;
    options.tensorOrder = 1;
    options.uTrans = 0;
    options.U = zonotope([0,0]);
    options.x0=center(zono0);
    options.R0=zono0;

    % compute reachable set
    [reach, rs] = timeSeries(tm0, @ll, options.timeStep, options.tFinal, optns);

    % simulation
    llSys = nonlinearSys(7,1,@(t, x, u) ll(x), options);
    %--------------------------------------------------------------------------
    runs = 60;
    fractionVertices = 0.5;
    fractionInputVertices = 0.5;
    inputChanges = 6;
    simRes = simulate_random(llSys, options, runs, fractionVertices, fractionInputVertices, inputChanges);

    %plot results--------------------------------------------------------------
    figure;
    hold on

    for i=1:rs
        %plot flowpipes
        zono = zono_of_taylm(reach{i}{2}, [names', {'t'}]);
        plotFilled(zono,options.projectedDimensions,[.5 .5 .5],'EdgeColor','black');
    end
    
    %plot initial set
    plotFilled(options.R0,options.projectedDimensions,'w','EdgeColor','black');

    %plot discrete sets
    for i=1:rs
        zono = zono_of_taylm(reach{i}{1}, names');
        plotFilled(zono,options.projectedDimensions,[.8 .8 .8],'EdgeColor','black');
    end
    
    %plot simulation results      
    for i=1:length(simRes.t)
        plot(simRes.x{i}(:,options.projectedDimensions(1)),simRes.x{i}(:,options.projectedDimensions(2)),'Color',0*[1 1 1]);
    end

end

function res = ll(x)
    res = [1.4*x(3) - 0.9*x(1);
        2.5*x(5)-1.5*x(2);
        0.6*x(7) - 0.8*x(2)*x(3);
        2 - 1.3*x(3)*x(4);
        0.7*x(1) - x(4)*x(5);
        0.3*x(1) - 3.1*x(6);
        1.8*x(6) - 1.5*x(2)*x(7)];
end