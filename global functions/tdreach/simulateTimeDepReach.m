
function [xs, ys] = simulateTimeDepReach(ivl, Nsamples, f, h, template, Nsteps)
    samples = sample_points(taylm(ivl), Nsamples);
    [~,Nsamples] = size(samples);
    simulations = zeros(2, Nsamples, Nsteps);
    for i = 1:Nsteps
        F_fitted = dependentTime(f, samples, h, template);
        for s = 1:Nsamples
            sample = samples(:,s);
            simulations(:,s,i) = sample;
            t = templatePoly(F_fitted, sample');
            [~, x, ~] = ode45(@(t, x) f(x), [0, t], sample);
            samples(:,s) = x(end,:);
        end
    end
    xs = simulations(1,:);
    ys = simulations(2,:);
end