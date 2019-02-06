function res = dependentTime(f, samples, h)
    % find closest point
    [n, N] = size(samples); % N sample points of dimension n
    
    f0xs = zeros(n, N);
    f1xs = zeros(n, N);
    x1s = zeros(n, N);
    for i = 1:N
        x0 = samples(:,i);
        f0xs(:,i) = f(x0);
        [~, x] = ode45(@(t, x) f(x), [0, h], x0);
        x1 = x(end,:)';
        x1s(:,i) = x1;
        f1xs(:,i) = f(x1);
    end
    nrml = mean(f0xs, 2);
    x0 = mean(samples(:,i),2);
    nrml = f(x0);
    c0s = zeros(1, N);
    c1s = zeros(1, N);
    for i = 1:N
        c0s(i) = nrml' * samples(:,i);
        c1s(i) = nrml' * x1s(:,i);
    end
    c = max(max(c0s(:)), min(c1s(:)));

    % find an optimal separating hyperplane (not a good idea)
    % (http://www.robots.ox.ac.uk/~az/lectures/ml/matlab2.pdf)
%    H = eye(n+1);
%    H(n+1,n+1) = 0;
%    A = -[eye(N), zeros(N); zeros(N), -eye(N)] * [[samples'; x1s'], ones(2*N,1)];
%    b = -ones(2*N,1);
%    opt = quadprog(H, zeros(n+1,1), A, b);
%    c = -opt(n + 1) - 1;
%   nrml = opt(1:n);
    
%   dists = zeros(N, 1);
%    for i = 1:N
%   x0 = samples(:,i);
%        [~, xs] = ode45(@(t, x) f(x), [0, h], x0);
%        x1 = xs(end,:)';
%        dists(i) = norm(x1 - x0);
%    end
%    [~, minidx] = min(dists); % Apparently a bad idea, it is a point that is furthest away from the limit behaviour.

    % assuming that all the others need shorter time to get there...
    % Might want to adjust the timestep a bit, though
%    [c, nrml] = normalHyperplaneAtStep(f, x0, h2); 
    
    explanatory = samples';
    response=zeros(N, 1);
    
    for i = 1:N
        r = returnTime(f, samples(:,i), c, nrml, 4*h);
        if isempty(r)
            disp(['Does not return within 100 times the original time,',...
                'something might be utterly wrong.'])
            response(i) = h;
        else
            response(i) = r;
        end
    end
%    Dr = approxReturnTimeDerivative(f, x0, h2);
    notimedep = [h, zeros(1, n)];
    res = nlinfit(explanatory,response,@templateLinear,notimedep);
    
    % check if the result makes sense, if not, smoothen the whole thing.
    ts = templateLinear(res, explanatory);
    tmin = min(ts);
    tmax = max(ts);
    while tmin <= 0 || tmax > 1
        disp ('Strange')
        res = (res + notimedep)/2;
    end
end
