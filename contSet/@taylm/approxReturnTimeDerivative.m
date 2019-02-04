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