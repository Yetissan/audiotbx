function [y] = autbx_stavgamp(x, window_type, window_length)
    w = window(window_type, window_length);
    y = conv(w, abs(x));

    y = y(1:length(x));
endfunction
