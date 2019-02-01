function [reach, rs] = timeSeries(x, f, h, T, optns)
    N = ceil(T / h);
    reach = cell(N, 2);
    t = 0.0;
    i = 1;
    while t < T
        [flowpipe, x] = timeStep(x, f, h, optns);
        reach{i}{1} = x;
        reach{i}{2} = flowpipe;
        t = t + h
        i = i + 1;
    end
    rs = i - 1;
end

function [flowpipe, final] = timeStep(x, f, h, optns)
    flowpipe = certify_step(f, x, h, optns);
    timestep = taylm(interval(1, 1), x(1).max_order);
    final = horner(flowpipe, {'t'}, timestep);
end
