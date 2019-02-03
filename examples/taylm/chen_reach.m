function chen_reach() % actually, now it is the example by
%   Neher et al: "On Taylor Model Based integration of ODEs"
    
    % optns for Picard iteration
    optns.picard_threshold = 1e-10;
    optns.picard_order = 6;
    optns.widening_scale = 1.1;
    optns.narrowing_scale = 1.1;
    optns.time_var = {'t'};
    optns.remainder_estimation = [0.0001; 0.0001];
    optns.parallelotope_factor = 1.1;

    % initalization of simulation and tdreach
    options.timeStep=0.3;
    options.tFinal=1.5;

    tm0 = taylm(interval([0.95; -1.05], [1.05; -0.95]),6, {'a'; 'b'}, 'int');
    zono0 = zono_of_taylm(tm0, ['a', 'b']);
    options.projectedDimensions=[1 2];

    % fixed options for simulation
    options.tStart=0;
    options.tensorOrder = 1;
    options.uTrans = 0;
    options.U = zonotope([0,0]);
    options.x0=center(zono0);
    options.R0=zono0;
    
    % function
    chen = @(x) [x(2); (x(1)^2)];
    
    % compute reachable set
    [reach, rs] = timeSeries(tm0, chen, options.timeStep, options.tFinal, optns);
    
    % simulation
    lorenz_sys = nonlinearSys(2,1,@chenEq, options);
    %--------------------------------------------------------------------------
    runs = 60;
    fractionVertices = 0.5;
    fractionInputVertices = 0.5;
    inputChanges = 1;
    simRes = simulate_random(lorenz_sys, options, runs, fractionVertices, fractionInputVertices, inputChanges);

    %plot results--------------------------------------------------------------
    plotOrder = 20;
    figure;
    hold on
    
    for i=1:rs
        %plot flowpipes
        zono = zono_of_taylm(reach{i}{2}, ['a', 'b', 't']);
        plotFilled(zono,options.projectedDimensions,[.5 .5 .5],'EdgeColor','black');
    end
    
    %plot initial set
    plotFilled(options.R0,options.projectedDimensions,'w','EdgeColor','black');

    %plot discrete sets
    for i=1:rs
        zono = zono_of_taylm(reach{i}{1}, ['a', 'b', 't']);
        plotFilled(zono,options.projectedDimensions,[.8 .8 .8],'EdgeColor','black');
    end
    
    %plot discrete as grid
    for i=1:rs
        i
        grid = grid_of_taylm(reach{i}{1}, options.projectedDimensions, 5);
        [sg, ~] = size(grid);
        for j = 1:sg
            zono = zono_of_taylm(grid{j}, ['a', 'b', 't']);
            plotFilled(zono,[1 2], [1.0 0.5 1.0],'EdgeColor','none');
        end
    end
    
    %plot simulation results      
    for i=1:length(simRes.t)
        plot(simRes.x{i}(:,options.projectedDimensions(1)),simRes.x{i}(:,options.projectedDimensions(2)),'Color',1*[1 0 1]);
    end
    
end

function [dx]=chenEq(~,x,~)
dx(1,1)=x(2);
dx(2,1)=x(1)^2;
end