function [y] = autbx_stenergy(x, window_type, window_length)
    w = window(window_type, window_length);
    y = conv(w.*w, x.*x);
    
    y = y(1:length(x));
endfunction
