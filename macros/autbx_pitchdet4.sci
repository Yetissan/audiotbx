function [y, t] = autbx_pitchdet4(x, n_start, n_end, fs, frame_size, stepping, fmin, fmax, num_of_freq_slots)
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
    
    if ((fmin < 0) | (fmax < 0) | (fmin >= fmax)) then
        error('Invalid frequency range.');
        return;
    end

    x = (x(:))';
    r = modulo(length(x) + frame_size - 1, frame_size);
    if (r ~= 0) then
        x = [x, zeros(1, frame_size - r)];
    end
    N = length(x);

    freq_slot_size = (fmax - fmin) ./ (num_of_freq_slots);
    freq_candidates = [];
    opt_freq_seq = [];
    prev_cost = [];
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
        
        tau_min = 1 ./ fmax;
        tau_max = 1 ./ fmin;
        M = max(cmndf);
        cmndf1 = cmndf;
        cmndf1(1 : max([1, floor(tau_min)])) = M;
        cmndf1(min([N, floor(tau_max)]) : $) = M;
        if (~isempty(find(cmndf1 < absthd))) then
            tmp = cmndf1;
            tmp(find(cmndf1 > absthd)) = absthd + 1;
            [m, idx1] = min(tmp);
        else
            [m, idx1] = min(cmndf1);
        end
        
        cmndf2 = cmndf;
        cmndf2(1 : max([1, floor(tau_min)])) = M;
        cmndf2(min([N, floor(0.75 .* tau_max)]) : $) = M;
        if (~isempty(find(cmndf2 < absthd))) then
            tmp = cmndf2;
            tmp(find(cmndf2 > absthd)) = absthd + 1;
            [m, idx2] = min(tmp);
        else
            [m, idx2] = min(cmndf2);
        end
        
        cmndf3 = cmndf;
        cmndf3(1 : max([1, floor(1.25 .* tau_min)])) = M;
        cmndf3(min([N, floor(tau_max)]) : $) = M;
        if (~isempty(find(cmndf3 < absthd))) then
            tmp = cmndf3;
            tmp(find(cmndf3 > absthd)) = absthd + 1;
            [m, idx3] = min(tmp);
        else
            [m, idx3] = min(cmndf3);
        end
        
        f1 = fs ./ idx1(1);
        f2 = fs ./ idx2(1);
        f3 = fs ./ idx3(1);
        if ((f1 > fmax) | (f1 < fmin)) then
            freq_candidates = [0, 0, 0; freq_candidates];
        else
            freq_candidates = [f1, f2, f3; freq_candidates];
        end
        
        if (i >= 2) then
            
        else
            prev_cost = 
        end
    end

    s = size(y);
    t = linspace((n_start - 1) / fs, (n_end - 1) / fs, s(2));
endfunction


