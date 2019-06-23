function [y] = autbx_threelevelclip(x, threshold)
    y = [];
    
    if (threshold <= 0) then
        error('Threshold must be greater than zero. ');
        return;
    end
    
    y = x;
    y(intersect(find(y >= -threshold), find(y <= threshold))) = 0;
    y(find(y > threshold)) = 1;
    y(find(y < -threshold)) = -1;
endfunction
