function res = dependentTime(f, samples, h, template)
    timeshorteningFactor = 0.8;
    % find closest point
    [n, N] = size(samples); % N sample points of dimension n
    dists = zeros(N, 1);
    for i = 1:N
        x0 = samples(:,i);
        [~, xs] = ode45(@(t, x) f(x), [0, h], x0);
        x1 = xs(end,:)';
        dists(i) = norm(x1 - x0);
    end
    [~, minidx] = min(dists);
    h2 = timeshorteningFactor * h;
    [c, nrml] = normalHyperplaneAtStep(f, samples(:,minidx), h2); % assuming that all the others need shorter time to get there... Might want to adjust the timestep a bit, though
    explanatory = samples';
    response=zeros(N, 1);
    for i = 1:N
        r = returnTime(f, samples(:,i), c, nrml, h);
        if isempty(r)
            error('Does not return within original time, try setting a lower timeshorteningFactor in dependentTime!')
        end
        response(i) = r;
    end
    Dr = approxReturnTimeDerivative(f, x0, h2);
    res = nlinfit(explanatory,response,template,[0, Dr, 0, 0, 0]);
end
