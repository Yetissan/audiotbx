funcprot(0);

function [y, t] = pitchdet(x, n_start, n_end, fs, frame_size, stepping, absthd, fmin, fmax)
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
    
    if (absthd <= 0) then
        error('absthd must be positive.');
        return;
    end
    
    if (stepping <= 0) then
        error('stepping must be positive.');
        return;
    end

    x = (x(:))';
    r = modulo(length(x) + frame_size - 1, frame_size);
    if (r ~= 0) then
        x = [x, zeros(1, frame_size - r)];
    end
    N = length(x);
    
    for i = 1 : stepping : (N - 2 .* frame_size + 2)
        s1 = x(i : i + frame_size - 1);

        // Average mean difference function
        amdf = [];
        for j = 0 : frame_size - 1
            s2 = x(i + j : i + frame_size + j - 1);
            amdf = [amdf, sum((s1 - s2) .^ 2)];
        end
        
        // Cumulative mean normalized difference function
        cmndf = [1];
        for k = 2 : length(amdf) / 2
            cmndf = [cmndf, amdf(k) .* (k - 1) ./ sum(amdf(2 : k))];
        end
        
        if (~isempty(find(cmndf < absthd))) then
            tmp = cmndf;
            tmp(find(cmndf > absthd)) = absthd + 1;
            [m, idx] = min(tmp);
        else
            [m, idx] = min(cmndf);
        end
        
        y = [y, fs/idx(1)];
    end

    s = size(y);
    t = linspace((n_start - 1) / fs, (n_end - 1) / fs, s(2));
endfunction

