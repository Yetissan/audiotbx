function [y, t, f] = autbx_spectogram2(x, fs, stft_size, window_type, ovlp_length)
    y = [];
    t = [];
    f = [];

    N = length(x);
    x = (x(:))';

    if (stft_size  > N) then
        error('stft_size must be smaller or equal to length of x.');
        return;
    end

    if ((ovlp_length > stft_size) | (ovlp_length < 0)) then
        error('Improper ovlp_length.');
        return;
    end
    
    r = modulo(N + ovlp_length, stft_size);
    if (r ~= 0) then
        x = [x, zeros(1, stft_size - r)];
        N = length(x);
    end

    w = window(window_type, stft_size);
    i = 1;
    while (i < N - stft_size)
        X = fft(x(i : min([i+stft_size-1, N])).*w);
        X = abs(X(2:stft_size/2)) ./ stft_size .* 2;
        y = [y; (X(:))'];

        i = i + (stft_size - ovlp_length);
    end

    s = size(y);
    t = linspace(0, N/fs, s(1));
    f = linspace(0, fs/2, s(2));
endfunction
