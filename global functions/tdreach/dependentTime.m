function res = dependentTime(f, samples, h)
    % find closest point
    [~, N] = size(samples); % N sample points of dimension ~
    dists = zeros(N, 1);
    
    x0 = mean(samples,2);

    
%    for i = 1:N
%   x0 = samples(:,i);
%        [~, xs] = ode45(@(t, x) f(x), [0, h], x0);
%        x1 = xs(end,:)';
%        dists(i) = norm(x1 - x0);
%    end
%    [~, minidx] = min(dists); % Apparently a bad idea, it is a point that is furthest away from the limit behaviour.
    
    
    h2 = 0.5 * h;
    % assuming that all the others need shorter time to get there...
    % Might want to adjust the timestep a bit, though
    [c, nrml] = normalHyperplaneAtStep(f, x0, h2); 
    explanatory = samples';
    response=zeros(N, 1);
    
    for i = 1:N
        r = returnTime(f, samples(:,i), c, nrml, 2*h);
        if isempty(r)
            error('Does not return within 100 times the original time, something is utterly wrong.')
        end
        response(i) = r;
    end
    Dr = approxReturnTimeDerivative(f, x0, h2);
    res = nlinfit(explanatory,response,@templateLinear,[h2, Dr]);
    
    % check if the result makes sense, if not, smoothen the whole thing.
    ts = templateLinear(res, explanatory);
end
