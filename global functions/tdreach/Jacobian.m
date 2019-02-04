function res = Jacobian(x, h, F)
    h = h*ones(size(x));
    res = (F(repmat(x,size(x'))+diag(h))-F(repmat(x,size(x'))))./h';
end
