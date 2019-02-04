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

function [position,isterminal,direction] = crossingEvent(c, n, x)
    position = n' * x - c; % The value that we want to be zero
    isterminal = 1;  % Halt integration 
    direction = 0;
end
