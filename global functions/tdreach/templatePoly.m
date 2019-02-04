function F = templatePoly(x, xdata)
    [n, ~] = size(xdata);
    F = xdata(:,1);
    for i = 1:n
        F(i, 1) = x(2).*xdata(i,1)   + x(3).*xdata(i,2) + ...
            x(4).*xdata(i,1)^2 + x(5).*xdata(i,1).*xdata(i,2) + x(6).*xdata(i,2)^2 + ...
            x(1);
    end
end
