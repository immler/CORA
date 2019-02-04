function [c, n] = normalHyperplaneAtStep(f, x0, h)
    [~, x] = ode45(@(t, x) f(x), [0, h], x0);
    x1 = x(end,:)';
    n = f(x1);
    c = f(x1)' * x1;
end

