function [y, t] = autbx_pitchdet4(x, n_start, n_end, fs, frame_size, stepping, fmin, fmax)
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

    freq_candidates = [];
    opt_freq_seq = [];
    prev_cost = ones(1, 3) .* %eps;
    seg_idx = 1;
    frame_idx = 1;
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

        // Search for the smallest delay which gives a local minimum of CMNDF
        // smaller than absthd in the following search ranges respectively:
        //
        // [tau_min, tau_max]           [fmin, fmax]            f1
        // [tau_min, 0.75 * taumax]     [1.33 * fmin, fmax]     f2
        // [1.25 * taumin, tau_max]     [fmin, 0.75 * fmax]     f3
        //
        // where fk (k = 1, 2, 3) be set to -1 if no eligible delay was found.
        // Especially if no f1 was found, we will take the frame being observed
        // by now as an unvoiced / unwanted frame thus setting f1, f2 and f3
        // to -1 as an indication of the fact.
        cmndf_lwr_envl = autbx_lwr_envelope(1 : length(cmndf), cmndf, 1);
        efctv_idx = find(cmndf_lwr_envl(2, :) < absthd);
        efctv_cmndf_lwr_envl_idx = cmndf_lwr_envl(1, efctv_idx);
        efctv_cmndf_lwr_envl_val = cmndf_lwr_envl(2, efctv_idx);
        if (length(efctv_cmndf_lwr_envl_idx) >= 2) then
            idx1 = efctv_cmndf_lwr_envl_idx(2);
            f1 = fs ./ idx1;
        else
            idx1 = -1;
            f1 = -1;
        end

        if (idx1 <> -1) then
            cmndf2 = cmndf;
            cmndf2(1 : min([length(cmndf2), max([1, floor(tau_min)])])) = M;
            cmndf2(max([1, min([length(cmndf2), floor(0.75 * tau_max)])]) : $) = M;
            cmndf2_lwr_envl = autbx_lwr_envelope(1 : length(cmndf2), cmndf2, 1);
            efctv_idx2 = find(cmndf2_lwr_envl(2, :) < absthd);
            efctv_cmndf2_lwr_envl_idx = cmndf2_lwr_envl(1, efctv_idx2);
            efctv_cmndf2_lwr_envl_val = cmndf2_lwr_envl(2, efctv_idx2);
            if (length(efctv_cmndf2_lwr_envl_idx) >= 2) then
                idx2 = efctv_cmndf2_lwr_envl_idx(2);
                f2 = fs ./ idx2;
            else
                idx2 = -1;
                f2 = -1;
            end

            cmndf3 = cmndf;
            cmndf3(1 : min([length(cmndf3), max([1, 1.25 * floor(tau_min)])])) = M;
            cmndf3(max([1, min([length(cmndf3), floor(tau_max)])]) : $) = M;
            cmndf3_lwr_envl = autbx_lwr_envelope(1 : length(cmndf3), cmndf3, 1);
            efctv_idx3 = find(cmndf3_lwr_envl(2, :) < absthd);
            efctv_cmndf3_lwr_envl_idx = cmndf3_lwr_envl(1, efctv_idx3);
            efctv_cmndf3_lwr_envl_val = cmndf3_lwr_envl(2, efctv_idx3);
            if (length(efctv_cmndf3_lwr_envl_idx) >= 2) then
                idx3 = efctv_cmndf3_lwr_envl_idx(2);
                f3 = fs ./ idx3;
            else
                idx3 = -1;
                f3 = -1;
            end
        else
            idx2 = -1;
            f2 = -1;
            idx3 = -1;
            f3 = -1;
        end
        
        state_cost = [f1, f2, f3] ./ M;

        if (i >= 2) then

        else
            prev_cost =
        end
        
        frame_idx = frame_idx + 1;
    end
    
    num_of_frames = frame_idx - 1;
    if (seg_idx < num_of_frames) then
    end

    s = size(y);
    t = linspace((n_start - 1) / fs, (n_end - 1) / fs, s(2));
endfunction


