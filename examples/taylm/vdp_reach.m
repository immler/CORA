function vdp_reach()
    % optns for Picard iteration
    optns.picard_threshold = 1e-5;
    optns.picard_order = 5;
    optns.widening_scale = 1.1;
    optns.narrowing_scale = 1.01;
    optns.time_var = {'t'};
    optns.remainder_estimation = [1e-5; 1e-5];
    optns.parallelotope_factor = 1.0;
    optns.shrinking_mod=1;
    optns.zonotope_enclosure_mod=20;
    optns.shrink_wrap_options = { 'small_factor' 1.01;
        'iter_max' 5;
        'q_max' 1.01;
        'q_tol' 1e-12};
    optns.timedependency = 1;
    optns.dt_mod = 1;
    
    % iniitalization of simulation and tdreach
    options.timeStep=0.1;
    options.tFinal=0.1;
    
    tm0 = taylm(interval([1.25; 2.25], [1.55; 2.35]),5, {'x'; 'y'}, 'int', 1e-12, 1e-8);
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
    
    [xs, ys] = simulateTimeDepReach(interval([1.25; 2.25], [1.55; 2.35]), 10, vdpt, 0.1, @templateLinear, 4);
    figure;
    plot(xs, ys, '.')
    return
    % experimenting with return times;;;
    if 1
        x0 = [1.4; 2.3];
        h = 0.1;
        dR = approxReturnTimeDerivative(vdpt, x0, h/2);
        xs = 1.25:0.05:1.55;
        ys = 2.25:0.01:2.35;
        [X, Y] = meshgrid(xs, ys);
        [n, m] = size(X);
        xp = zeros(n,m);
        yp = zeros(n,m);
        f = vdpt;
        samples = [X(:), Y(:)]';
        [~,Nsamples] = size(samples);
        Nsteps = 70;
        simulations = zeros(2, Nsamples, Nsteps);
        for i = 1:Nsteps
            F_fitted = dependentTime(f, samples, h, @templateLinear);
            for s = 1:Nsamples
                sample = samples(:,s);
                simulations(:,s,i) = sample;
                t = templatePoly(F_fitted, sample');
                [~, x, ~] = ode45(@(t, x) vdpt(x), [0, t], sample);
                samples(:,s) = x(end,:);
            end
        end
        % x: coefficients of polynomial encoding time dependency
        % xdata: input values
        % ydata: Poincare map onto the thing
        xs = zeros(Nsteps*Nsamples,1);
        ys = zeros(Nsteps*Nsamples,1);
        for n=1:Nsteps
            for s=1:Nsamples
                for i=1:2
                    
                end
            end
        end
        figure
        plot(simulations(1,:), simulations(2,:), 'o')
        return%:)
    end

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
%        plot(simRes.x{i}(:,options.projectedDimensions(1)),simRes.x{i}(:,options.projectedDimensions(2)),'Color',0*[1 1 1]);
    end

    % plot originial fitting
    plot(xf, yf, 'o')

end

