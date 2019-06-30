function [y, t] = autbx_pitchdet1(x, n_start, n_end, fs, window_type, frame_size, ovlp_length, medfilt_order)
    y = [];
    t = [];

    if ((n_start <= 0) | (n_end > length(x))) then
        error('Out of bound.');
        return;
    end

    x = x(n_start:n_end);
    N = length(x);

    if (frame_size  > N) then
        error('frame_size must be smaller or equal to length of x.');
        return;
    end

    if ((ovlp_length > frame_size) | (ovlp_length < 0)) then
        error('Improper ovlp_length.');
        return;
    end

    x = (x(:))';
    r = modulo(length(x) + ovlp_length, frame_size);
    if (r ~= 0) then
        x = [x, zeros(1, frame_size - r)];
    end
    N = length(x);
    
    w = window(window_type, frame_size);
    wcr = ifft(abs(fft(w)) .^ 2);
    wcr = wcr ./ wcr(1);
    
    i = 1;
    while (i < N - frame_size)
        s = x(i : min([i + frame_size - 1, N]));
        s = s - mean(s);
        
        s = (s(:))';
        ss = s .* w;
        cr = corr(ss, ss, frame_size) ./ wcr;
        
        // Find first local minimum in cr.
        j = 1;
        while (j + 2 < frame_size / 2)
            if ((cr(j) <= cr(j+1)) & (cr(j+1) <= cr(j+2))) then
                break;
            end

            j = j + 1;
        end
        [p, k] = max(cr(j+1:$));
        y = [y, fs / (k+j)];
            
        i = i + (frame_size - ovlp_length);
    end
    
    y = autbx_medfilter(y, medfilt_order);

    z = size(y);
    t = linspace((n_start - 1)/fs, (n_end - 1)/fs, z(2));
endfunction

