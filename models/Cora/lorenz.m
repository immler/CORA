function [dx] = lorenz(~, x, ~)

beta=8.0/3.0;
rho=28;
sigma=10;

dx(1,1) = sigma * (x(2) - x(1));
dx(2,1) = x(1) * (rho - x(3)) - x(2);
dx(3,1) = x(1) * x(2) - beta * x(3);

end

