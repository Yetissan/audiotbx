function [y] = autbx_stzerocross(x, window_type, window_length)
    w = window(window_type, window_length) ./ (2*window_length);
    
    x = (x(:))';
    n = length(x);
    
    x1 = [];
    for i = 1:n
        if (x(i) >= 0) then
            x1(i) = 1;
        else
            x1(i) = -1;
        end
    end
    
    x1 = (x(:))';
    x2 = [0, x1(1:$-1)];
    
    y = conv(abs(x1-x2), w);
    y = y(1:length(x));
endfunction

