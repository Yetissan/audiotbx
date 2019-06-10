function [X] = autbx_cepstrum(x, window_type)
    X = ifft(log(fft(x .* window(window_type, length(x)))));
endfunction

