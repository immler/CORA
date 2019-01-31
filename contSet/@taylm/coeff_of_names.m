function res = coeff_of_names ( obj, names, ps )
    indices = arrayfun (@(n) index_of(obj, n), names);
    deg = sum(ps);
    if all(arrayfun(@is_pos, indices))
        [mons, ~] = size(obj.monomials);
        for m = 1:mons
            if obj.monomials(m, 1) == deg
                mismatch = 0;
                for n = 1:size(names)
                    if obj.monomials(m, indices(n)) ~= ps(n)
                        mismatch = 1;
                    end
                end
                if mismatch == 0
                    res = obj.coefficients(m);
                    return
                end
            end
        end
        res = 0;
    else
        res = 0;
    end
end

function res = is_pos(n)
  if n > 0
      res = 1;
  else
      res = 0;
  end
end