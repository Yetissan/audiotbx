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

    if ((fmin < 0) | (fmax < 0) | (fmin >= fmax) | (fmax >= fs ./ 2)) then
        error('Invalid range of frequency.');
        return;
    end

    x = (x(:))';
    r = modulo(length(x) + frame_size - 1, frame_size);
    if (r ~= 0) then
        x = [x, zeros(1, frame_size - r)];
    end
    N = length(x);

    freq_candidates = [];
    opt_freq_idx = [];
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

        // Search for the smallest delay which gives a minimum of CMNDF
        // smaller than absthd in the following search ranges respectively:
        //
        // [tau_min, tau_max]           [fmin, fmax]            f1
        // [tau_min, 0.75 * taumax]     [1.33 * fmin, fmax]     f2
        // [1.25 * taumin, tau_max]     [fmin, 0.75 * fmax]     f3
        //
        // where fk (k = 1, 2, 3) be set to 0 if no eligible delay was found.
        // Especially if no feasible f1 was found, we will take the frame being
        // observed by now as an unwanted frame (happens to be an unvoiced one
        // sometimes) thus setting f1, f2 and f3 to 0 as an indication of the fact.
        tau_min = fs ./ fmax;
        tau_max = fs ./ fmin;
        cmndf_max = max(cmndf);

        search_ranges = [   tau_min,        tau_max;            ...
                            tau_min,        0.75 * tau_max;     ...
                            1.25 * tau_min, tau_max ];
        search_ranges_size = size(search_ranges);
        num_of_freq_states = search_ranges_size(1);
        curr_freq_cndds = zeros(1, num_of_freq_states);
        curr_delay_cndds = zeros(1, num_of_freq_states);
        for j = 1 : num_of_freq_states
            curr_search_range = search_ranges(i, :);
            __cmndf = cmndf;
            __cmndf(find(__cmndf >= absthd)) = cmndf_max;
            [m, __idx] = min(__cmndf);
            __idx = __idx(1);
            if ((__idx >= curr_search_range(1)) & (__idx <= curr_search_range(2))) then
                curr_freq_cndds(j) = fs ./ __idx;
                curr_delay_cndds(j) = __idx;
            elseif (j == 1) then
                break;
            end
        end
        freq_candidates = [freq_candidates; curr_freq_cndds];

        if (frame_idx == 1) then
            prev_total_full_cost = ones(1, num_of_freq_states);
            prev_avail_freqs_idx = find(curr_freq_cndds > 0);
            prev_avail_freqs_idx = (prev_avail_freqs_idx(:))';
            prev_freq_cndds = curr_freq_cndds;
        end

        if (frame_idx >= 2) then
            curr_avail_freqs_idx = find(curr_freq_cndds > 0);
            if ((~isempty(prev_avail_freqs_idx)) & (~isempty(curr_avail_freqs_idx))) then
                // $$ J_{T}(i;j,k) = \left| \frac{1}{\ln{2 f_{i-1}^{j}}} \cdot \left( \ln{f_{i}^{k}} - \ln{f_{i-1}^{j}}  \right)\right|  $$
                transition_costs = zeros(num_of_freq_states, num_of_freq_states);
                for p = find(prev_freq_cndds > 0)
                    for q = find(curr_freq_cndds > 0)
                        transition_costs(p, q) = abs((log(curr_freq_cndds(q)) - log(prev_freq_cndds(p))) ./ log(2 * prev_freq_cndds(p)));
                    end
                end

                // $$ J_{S}(i;k) = \frac{d^{\prime}(\tau_{k})}{\max d^{\prime}(\cdot)} $$
                state_costs = -1 * ones(1, num_of_freq_states);
                for p = find(curr_delay_cndds > 0)
                    state_costs(p) = cmndf(curr_delay_cndds) ./ cmndf_max;
                end

                // $$ J_{F}(i;j,k) = \alpha J_{T}(i;j,k) + (1 - \alpha) J_{S}(i;k)  $$
                alpha = 0.5;
                full_costs = -1 * ones(num_of_freq_states, num_of_freq_states);
                for p = find(prev_freq_cndds > 0)
                    for q = find(curr_freq_cndds > 0)
                        full_costs(p, q) = (alpha .* transition_costs(p, q)) + ((1 - alpha) .* state_costs(q));
                    end
                end

                normalize_factor = max(full_costs);
                if (normalize_factor <> -1) then
                    for p = find(prev_freq_cndds > 0)
                        for q = find(curr_freq_cndds > 0)
                            full_costs(p, q) = full_costs(p, q) ./ normalize_factor;
                        end
                    end
                end

                curr_total_full_cost = zeros(1, num_of_freq_states);
                for j = curr_avail_freqs_idx
                    cost0 = prev_total_full_cost .* full_costs(:, j)';
                    curr_opt_paths = intersect(prev_avail_freqs_idx, find(full_costs(:, j) > 0));

                    min_cost_path = curr_opt_paths(1);
                    if (length(curr_opt_paths) >= 2) then
                        for (k = curr_opt_paths(2 : $))
                            if (cost0(k) < cost0(min_cost_path)) then
                                min_cost_path = k;
                            end
                        end
                    end
                    curr_total_full_cost(j) = cost0(min_cost_path);
                    opt_freq_idx(i, j) = min_cost_path;
                end

                normalize_factor = max(curr_total_full_cost);
                curr_total_full_cost = curr_total_full_cost ./ normalize_factor;

                prev_avail_freqs_idx = curr_avail_freqs_idx;
                prev_total_full_cost = curr_total_full_cost;
                prev_freq_cndds = curr_freq_cndds;
            else
                prev_total_full_cost = ones(1, num_of_freq_states);
                prev_avail_freqs_idx = find(freq_candidates(frame_idx, :) > 0);
                prev_avail_freqs_idx = (prev_avail_freqs_idx(:))';
                prev_freq_cndds = curr_freq_cndds;
                opt_freq_idx = [opt_freq_idx; zeros(1, num_of_freq_states)];
            end
        end

        frame_idx = frame_idx + 1;
    end

    num_of_frames = frame_idx - 1;

    term_opt_freq_idx = 1;
    if (~isempty(prev_avail_freqs_idx)) then
        term_opt_freq_idx = prev_avail_freqs_idx(1);
        if (length(prev_avail_freqs_idx) >= 2) then
            for k = prev_avail_freqs_idx(2 : $)
                if (prev_total_full_cost(k) < prev_total_full_cost(term_opt_freq_idx)) then
                    term_opt_freq_idx = k;
                end
            end
        end
    end

    y = freq_candidates(num_of_frames, term_opt_freq_idx);
    prev_opt_freq_idx = term_opt_freq_idx;
    for k = (num_of_frames - 1) : -1 : 1
        curr_opt_freq_idx = opt_freq_idx(k + 1, prev_opt_freq_idx);
        if (curr_opt_freq_idx > 0) then
            y = [freq_candidates(k, curr_opt_freq_idx), y];
        else
            y = [0, y];
        end
        prev_opt_freq_idx = curr_opt_freq_idx;
    end

    s = size(y);
    t = linspace((n_start - 1) / fs, (n_end - 1) / fs, s(2));
endfunction


