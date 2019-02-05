function F = templateLinear(x, xdata)
    [n, ~] = size(xdata);
    F = xdata(:,1);
    for i = 1:n
        for d = 1:n
            F(i, 1) = x(1) + xdata(i,:) * x(2:n+1);
        end
    end
end
