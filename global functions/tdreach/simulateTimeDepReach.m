
function [ts, xs, ws] = simulateTimeDepReach(ivl, Nsamples, f, h, template, T, dodep)
    samples = sample_points(taylm(ivl), Nsamples);
    [D,Nsamples] = size(samples);
    ts = zeros(Nsamples,1);
    ws = zeros(Nsamples,1);
    t = 0;
    i = 0;
    while t < T
        i = i + 1;
        disp(['step no ', num2str(i), ' t = ', nu2str(t)])
        if dodep
            F_fitted = dependentTime(f, samples, h, template);
        end
        tmin = inf;
        for s = 1:Nsamples
            sample = samples(:,s);
            simulations(:,s,i) = sample;
            if dodep
                t = template(F_fitted, sample');
            else
                t = h;
            end
            tmin = min(t, tmin);
            [~, x, ~] = ode45(@(t, x) f(x), [0, t], sample);
            samples(:,s) = x(end,:);
        end
        t = t + tmin;
        ts(i) = t;
        % record maximum values
        wd = zeros(D, 1);
        for d = 1:D
            wd(d) = max(samples(d,:)) - min(samples(d,:));
        end
        ws(i) = max(wd);
    end
    for i = 1:D
        xs{i} = simulations(i,:);
    end
end