function res = remainders( obj )
    res = arrayfun (@(o) o.remainder, obj);
end
