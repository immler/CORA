function [ts, xs, ws] = simulateTimeDepReach(ivl, sample_gridN, f, h0, T, dodep)
    samples = sample_points(taylm(ivl), sample_gridN);
    [D,Nsamples] = size(samples);
    ts = zeros(Nsamples,1);
    ws = zeros(Nsamples,1);
    t = 0;
    i = 0;
    while t < T
        i = i + 1;
        disp(['step no ', num2str(i), ', t = ', num2str(t)])
        if dodep
            timedep = dependentTime(f, samples, h0);
        end
        timedeps = timedep' * [ones(1, Nsamples); samples];
        hmin = inf;
        for s = 1:Nsamples
            sample = samples(:,s);
            simulations(:,s,i) = sample;
            if dodep
                h = timedep' * [1; sample];
            else
                h = h0/4;
            end
            hmin = min(h, hmin);
            [~, x, ~] = ode45(@(t, x) f(x), [0, h], sample);
            samples(:,s) = x(end,:);
        end
        t = t + hmin;
        ts(i) = t;
        % record maximum values
        wd = zeros(D, 1);
        lower=zeros(D,1);
        upper=zeros(D,1);
        for d = 1:D
            lower(d) = min(samples(d,:));
            upper(d) = max(samples(d,:));
            wd(d) = max(samples(d,:)) - min(samples(d,:));
        end
        ws(i) = max(wd);
    end
    for i = 1:D
        xs{i} = simulations(i,:);
    end
end