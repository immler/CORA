function chen_reach() % actually, now it is the example by
%   Neher et al: "On Taylor Model Based integration of ODEs"
    
    % optns for Picard iteration
    optns.picard_threshold = 1e-5; % when are we happy with a picard approximation
    optns.picard_order = 3;
    optns.widening_scale = 1.1;
    optns.narrowing_scale = 1.1;
    optns.time_var = {'t'};
    optns.remainder_estimation = [0.001; 0.001];
    
    % iniitalization of simulation and tdreach
    options.tFinal=0.02;

    tm0 = taylm(interval([0.95; -1.05], [1.05; -0.95]),3, {'a'; 'b'}, 'int');
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
    certify_step(chen, tm0, 0.1, optns)
    
    return
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
    
end
