function [reach, rs] = timeSeries(x, f, h, T, optns)
    N = ceil(T / h);
    reach = cell(N, 2);
    t = 0.0;
    i = 0;
    while t < T
        [flowpipe, x] = timeStep(x, f, h, optns, i);
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
                x
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
        i = i + 1;
        reach{i}{1} = x;
        reach{i}{2} = flowpipe;
    end
    rs = i;
end

function res = approxReturnTimeDerivative(f, x0, h)
    JacobianDelta = 1e-2;
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

function [flowpipe, final] = timeStep(x, f, h, optns, i)
    flowpipe = certify_step(f, x, h, optns);
    % TODO: ensure timestep is included in [0, h]
    if optns.timedependency == 1 && mod(i, optns.dt_mod) == 0
        c = arrayfun(@getCoef, center(x));
        h2 = 0.5 * h;
        Dr = approxReturnTimeDerivative(f, c, h2)
%            timestep = 1 + 2 * (x - c) * Dr / h;
% s on [-1; 1]
% t on [0; h]
% t = (s + 1)*h/2
% s = 2*t/h - 1
        xwr = set_remainders(x, repmat(interval(0,0), size(x)))
        timestep = Dr * (xwr - c) / h % checked, this is perfectly right!
        itstp = interval(timestep)
        if interval(timestep) <= interval(-1, 1)
        else
            error ('The time is too wrong')
        end
    else
        timestep = taylm(interval(1, 1), x(1).max_order);
    end
    
    final = horner(flowpipe, {'t'}, timestep );
end
