function res = approxReturnTimeDerivative(f, x0, h)
    JacobianDelta = 1e-5;
    JacobianTimeFactor = 5;
    [~, m] = size(x0);
    if m ~= 1
        error('approx returntimederivative expecting column vector')
    end
    % make an approximate time step
    [c, n] = normalHyperplaneAtStep(f, x0, h);
    res = Jacobian(x0, JacobianDelta, @(x) returnTime(f, x, c, n, JacobianTimeFactor*5));
end