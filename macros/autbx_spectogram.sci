function [y, t, f] = autbx_spectogram(x, fs, stft_size, window_type)
    [y, t, f] = autbx_spectogram2(x, fs, stft_size, window_type, stft_size/2);
endfunction
