function [reach, rs] = timeSeries(x, f, h, T, optns)
    N = T / h + 1;
    reach = cell(N, 2);
    t = 0.0
    i = 1
    while t < T
        dh = min(h, T - t);
        [flowpipe, x] = timeStep(x, f, dh, optns);
        reach{i}{1} = x;
        reach{i}{2} = flowpipe;
        t = t + dh
        i = i + 1;
    end
    rs = i - 1
end

function [flowpipe, final] = timeStep(x, f, h, optns)
    flowpipe = picardIter(x, f, h, optns);
    timestep = taylm(interval(h, h));
    final = horner(flowpipe, 't', timestep);
end
