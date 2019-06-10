function [y] = autbx_centerclip(x, threshold)
    y = [];
    
    if (threshold <= 0) then
        error('Threshold must be greater than zero. ');
        return;
    end
    
    for i = 1:length(x)
        if (x(i) > threshold) then
            y(i) = x(i) - threshold;
        elseif (x(i) < -(threshold)) then
            y(i) = x(i) + threshold;
        else
            y(i) = 0;
        end
    end
endfunction
