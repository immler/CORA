function [reach, rs] = timeSeries(x, f, h, T, optns)
    N = ceil(T / h);
    reach = cell(N, 2);
    t = 0.0;
    i = 1;
    while t < T
        [flowpipe, x] = timeStep(x, f, h, optns);
        t = t + h
        w = max(arrayfun(@(x) rad(x.remainder), x))
        
        if mod(i, optns.zonotope_enclosure_mod) == 0
            x = order_1(x);
        end

        if mod(i, optns.shrinking_mod) == 0
            try
                tic
                x = shrink_wrap(x, optns.shrink_wrap_options);
                disp(['Shrink Wrap:', num2str(toc), ' s'])
            catch M
                disp('Catch')
                disp(M)
                tic
                wrap = parallelotope_wrap(x);
                disp(['Parallelotope Wrap:', num2str(toc), ' s'])
                if interval(wrap) <= optns.parallelotope_factor * interval(x)
                    x = wrap;
                elseif min(rad(interval(x))) >= 1
                    disp('Man, they are growing so fast')
                    break
                else
                    disp('No Parallelotope!')
                end
            end
        end
        reach{i}{1} = x;
        reach{i}{2} = flowpipe;
        i = i + 1;
    end
    rs = i - 1;
end

function [flowpipe, final] = timeStep(x, f, h, optns)
    flowpipe = certify_step(f, x, h, optns);
    
    switch optns.timedependency
        case 0
            timestep = taylm(interval(1, 1), x(1).max_order);    
        case 1
            c = arrayfun(@getCoef, center(x));
            Dr = approxReturnTimeDerivative(f, c, h);
            timestep = 1 + 2 * (x - c) * Dr / h;
    end
    
    final = horner(flowpipe, {'t'}, timestep );
end

function [position,isterminal,direction] = crossingEvent(c, n, x)
    position = n' * x - c; % The value that we want to be zero
    isterminal = 1;  % Halt integration 
    direction = 0;
end

function res = Jacobian(x, h, F)
    h = h*ones(size(x));
    res = (F(repmat(x,size(x'))+diag(h))-F(repmat(x,size(x'))))./h';
end

function te = returnTime(f, x0, c, n, T)
    [~, m] = size(x0);
    te = zeros(1,m);
    for i = 1:m
        te(1,i) = s_returnTime(f, x0(:,i), c, n, T);
    end
end

function te = s_returnTime(f, x0, c, n, T)
    Opt = odeset('Events', @(~, x) crossingEvent(c, n, x));
    [~, ~, te, ~, ~] = ode45(@(t, x) f(x), [0, T], x0, Opt);
end

function res = approxReturnTimeDerivative(f, x0, h)
    JacobianDelta = 1e-5;
    JacobianTimeFactor = 5;
    [~, m] = size(x0);
    if m ~= 1
        error('approx returntimederivative expecting column vector')
    end
    % make an approximate time step
    [~, x] = ode45(@(t, x) f(x), [0, h], x0);
    x1 = x(end,:)';
    n = f(x1);
    c = f(x1)' * x1;    
    res = Jacobian(x0, JacobianDelta, @(x) returnTime(f, x, c, n, JacobianTimeFactor*5));
end
