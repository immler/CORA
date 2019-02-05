function F = templateLinear(x, xdata)
    [n, D] = size(xdata);
    F = xdata(:,1);
    for i = 1:n
        for d = 1:D
            F(i, 1) = x(1) + x(2:D+1) * xdata(i,:)';
        end
    end
end
