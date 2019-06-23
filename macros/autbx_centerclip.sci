function [y] = autbx_centerclip(x, threshold)
    y = [];
    
    if (threshold <= 0) then
        error('Threshold must be greater than zero. ');
        return;
    end
    
    y = x;

    untouched = intersect(find(x >= -threshold), find(x <= threshold));
    greater_than_thrd = find(x > threshold);
    smaller_than_thrd = find(x < -threshold);
    y(greater_than_thrd) = y(greater_than_thrd) - threshold;
    y(smaller_than_thrd) = y(smaller_than_thrd) + threshold;
    y(untouched) = x(untouched);
endfunction
