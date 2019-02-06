function res = dependentTime(f, samples, h)
    % find closest point
    [n, N] = size(samples); % N sample points of dimension n
    
    f0xs = zeros(n, N);
    f1xs = zeros(n, N);
    n0xs = zeros(n, N);
    n1xs = zeros(n, N);
    x1s = zeros(n, N);
    for i = 1:N
        x0 = samples(:,i);
        fx0 = f(x0);
        f0xs(:,i) = fx0;
        n0xs(:,i) = fx0 / norm(fx0);
        [~, x] = ode45(@(t, x) f(x), [0, h], x0);
        x1 = x(end,:)';
        x1s(:,i) = x1;
        fx1 = f(x1);
        f1xs(:,i) = fx1;
        n1xs(:,i) = fx1 / norm(fx1);
    end
    lower0 = zeros(N,1);
    upper0 = zeros(N,1);
    lower1 = zeros(N,1);
    upper1 = zeros(N,1);
    for i = 1:N
        lower0(i) = max(samples' * n0xs(:,i));
        upper0(i) = min(x1s' * n0xs(:,i));
        lower1(i) = max(samples' * n1xs(:,i));
        upper1(i) = min(x1s' * n1xs(:,i));
    end
    [sep0, isep0] = max(upper0 - lower0);
    [sep1, isep1] = max(upper1 - lower1);
    sep = max(sep0, sep1);
    if max(sep0, sep1) <= 0
        disp('no separation for directions of flow -> trying to find a separating hyperplane!')
        % find an optimal separating hyperplane (not a good idea)
        % (http://www.robots.ox.ac.uk/~az/lectures/ml/matlab2.pdf)
        H = eye(n+1);
        H(n+1,n+1) = 0;
        A = -[eye(N), zeros(N); zeros(N), -eye(N)] * [[samples'; x1s'], ones(2*N,1)];
        b = -ones(2*N,1);
        opt = quadprog(H, zeros(n+1,1), A, b);
        c = -opt(n + 1) - 1;
        nrml = opt(1:n);
    elseif sep0 >= sep1
        nrml = n0xs(:,isep0);
        c = upper0(isep0);
    else
        nrml = n1xs(:,isep1);
        c = upper1(isep1);
    end
    
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
