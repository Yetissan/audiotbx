function [x] = autbx_icepstrum(X, window_type)
    x = ifft(exp(fft(X))) ./ window(window_type, length(X));
endfunction

