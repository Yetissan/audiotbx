function [y, t] = autbx_pitchdet4(x, n_start, n_end, fs, frame_size, ovlp_length, fmin, fmax, absthd)
    if ((n_start <= 0) | (n_end > length(x))) then
        error('Out of bound.');
        return;
    end

    x = x(n_start:n_end);
    N = length(x);

    if (frame_size  > N) then
        error('frame_size must be an integer smaller or equal to length of x.');
        return;
    end

    if ((ovlp_length < 0) | (ovlp_length > frame_size)) then
        error('ovlp_length must be a non-negative integer smaller than frame_size.');
        return;
    end

    if (absthd < 0) then
        error('absthd: A non-negative real number is required.');
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
    opt_freqs_idx = [];
    frame_idx = 1;
    last_valid_frm_idx = 0;
    last_avail_freqs_idx = [];
    last_avl_tot_full_cost = [];
    is_prev_frm_valid = 0;
    is_curr_frm_valid = 0;
    part_opt_freqs_idx = [];
    part_term_frms_idx = [];
    for i = 1 : (frame_size - ovlp_length) : (N - 2 .* frame_size + 2)
        printf('Frame #%d\t', frame_idx);
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
        // [tau_min, tau_max]           [fmin, fmax]                --> f1
        // [0.707f1, 1.414f1]           [1.414tau1, 2.828tau1]      --> f2
        // [1.414f1, 2.828f1]           [0.353tau1, 0.707tau1]      --> f3
        //
        // where fk (k = 1, 2, 3) be set to 0 if no eligible delay was found.
        // Especially if no feasible f1 was found, we will take the frame being
        // observed by now as an unwanted frame (happens to be an unvoiced one
        // sometimes) thus setting f1, f2 and f3 to 0 as an indication of the fact.
        tau_min = fs ./ fmax;
        tau_max = fs ./ fmin;
        cmndf_max = max(cmndf);

        search_ranges = [   tau_min,        tau_max;            ...
                            0,              0;                  ...
                            0,              0   ];
        search_ranges_size = size(search_ranges);
        num_of_freq_states = search_ranges_size(1);
        curr_freq_cndds = zeros(1, num_of_freq_states);
        curr_delay_cndds = zeros(1, num_of_freq_states);
        __pkreltol1 = 0.5;
        __pkreltol2 = 0.35;
        for j = 1 : num_of_freq_states
            curr_search_range = floor(search_ranges(j, :));
            __cmndf = cmndf;
            __cmndf(1 : min([length(__cmndf), max([1, floor(curr_search_range(1))])])) = cmndf_max;
            __cmndf(max([1, min([length(__cmndf), floor(curr_search_range(2))])]) : $) = cmndf_max;
            __cmndf(find(__cmndf >= (1 + __pkreltol1)*absthd)) = cmndf_max;
            __cmndf_lwr_envlp = autbx_lwr_envelope(1 : length(__cmndf), __cmndf, 1, 0);
            __cmndf_lwr_envlp_idx = __cmndf_lwr_envlp(1, :);
            __cmndf_lwr_envlp_val = __cmndf_lwr_envlp(2, :);
            if (~isempty(__cmndf_lwr_envlp_idx)) then
                __cmndf_lwr_envlp_idx = __cmndf_lwr_envlp_idx(find((__cmndf_lwr_envlp_idx >= curr_search_range(1)) & ...
                    (__cmndf_lwr_envlp_idx <= curr_search_range(2))));
                __bidx = __cmndf_lwr_envlp_idx(find(__cmndf_lwr_envlp_val <= (absthd * (1 + __pkreltol2))));
                if (~isempty(__bidx) & (__bidx(1) >= max([tau_min, curr_search_range(1)])) & ...
                    (__bidx(1) <= min([tau_max, curr_search_range(2)]))) then
                    __idx = min(__bidx);
                    curr_freq_cndds(j) = fs ./ __idx;
                    curr_delay_cndds(j) = __idx;
                    if (j == 1) then
                        is_curr_frm_valid = 1;
                        last_valid_frm_idx = frame_idx;
                        search_ranges(2, 1) = 1.414 * __idx;
                        search_ranges(2, 2) = 2.828 * __idx;
                        search_ranges(3, 1) = 0.353 * __idx;
                        search_ranges(3, 2) = 0.707 * __idx;
                    end
                else
                    if (j == 1) then
                        is_curr_frm_valid = 0;
                        break;
                    end
                end
            elseif (j == 1) then
                is_curr_frm_valid = 0;
                break;
            end
        end
        freq_candidates = [freq_candidates; curr_freq_cndds];

        printf('[ %d %d ]\n', is_prev_frm_valid, is_curr_frm_valid);
        select ([is_prev_frm_valid, is_curr_frm_valid])
            case [0, 0] then
                opt_freqs_idx = [opt_freqs_idx; zeros(1, num_of_freq_states)];
                is_prev_frm_valid = 0;
            case [0, 1] then
                curr_avail_freqs_idx = find(curr_freq_cndds > 0);
                curr_state_costs = ones(1, num_of_freq_states);
                for p = curr_avail_freqs_idx
                    curr_state_costs(p) = cmndf(curr_delay_cndds(p)) ./ cmndf_max + %eps;
                end
                normalize_factor = max(curr_state_costs);
                for p = curr_avail_freqs_idx
                    curr_state_costs(p) = curr_state_costs(p) ./ normalize_factor + %eps;
                end
                prev_total_full_cost = curr_state_costs;
                prev_avail_freqs_idx = curr_avail_freqs_idx;
                prev_freq_cndds = curr_freq_cndds;
                opt_freqs_idx = [opt_freqs_idx; zeros(1, num_of_freq_states)];
                last_valid_frm_idx = frame_idx;
                last_avail_freqs_idx = prev_avail_freqs_idx;
                last_avl_tot_full_cost = curr_state_costs;
                is_prev_frm_valid = 1;
            case [1, 0] then
                term_opt_freq_idx = prev_avail_freqs_idx(1);
                if (length(prev_avail_freqs_idx) >= 2) then
                    for k = prev_avail_freqs_idx(2 : $)
                        if (prev_total_full_cost(k) < prev_total_full_cost(term_opt_freq_idx)) then
                            term_opt_freq_idx = k;
                        end
                    end
                end
                part_opt_freqs_idx = [term_opt_freq_idx, part_opt_freqs_idx];
                part_term_frms_idx = [frame_idx, part_term_frms_idx];
                opt_freqs_idx = [opt_freqs_idx; zeros(1, num_of_freq_states)];
                is_prev_frm_valid = 0;
            case [1, 1] then
                curr_avail_freqs_idx = find(curr_freq_cndds > 0);
                // $$ J_{T}(i;j,k) = \left| \frac{1}{\ln{2 f_{i-1}^{j}}} \cdot \left( \ln{f_{i}^{k}} - \ln{f_{i-1}^{j}}  \right)\right|  $$
                transition_costs = ones(num_of_freq_states, num_of_freq_states);
                for p = prev_avail_freqs_idx
                    for q = curr_avail_freqs_idx
                        transition_costs(p, q) = ...
                            abs((log(curr_freq_cndds(q)) - log(prev_freq_cndds(p))) ./ log(2 * prev_freq_cndds(p))) + %eps;
                    end
                end

                // $$ J_{S}(i;k) = \frac{d^{\prime}(\tau_{k})}{\max d^{\prime}(\cdot)} $$
                state_costs = ones(1, num_of_freq_states);
                for p = curr_avail_freqs_idx
                    state_costs(p) = cmndf(curr_delay_cndds(p)) ./ cmndf_max + %eps;
                end

                // $$ J_{F}(i;j,k) = \alpha J_{T}(i;j,k) + (1 - \alpha) J_{S}(i;k)  $$
                alpha = 0.5;
                full_costs = -1 * ones(num_of_freq_states, num_of_freq_states);
                for p = prev_avail_freqs_idx
                    for q = curr_avail_freqs_idx
                        full_costs(p, q) = (alpha .* transition_costs(p, q)) + ((1 - alpha) .* state_costs(q)) + %eps;
                    end
                end

                normalize_factor = max(full_costs);
                if ((normalize_factor <> -1) & (normalize_factor <> 0)) then
                    for p = prev_avail_freqs_idx
                        for q = curr_avail_freqs_idx
                            full_costs(p, q) = full_costs(p, q) ./ normalize_factor;
                        end
                    end
                end

                curr_total_full_cost = zeros(1, num_of_freq_states);
                for j = curr_avail_freqs_idx
                    cost0 = prev_total_full_cost .* full_costs(:, j)';
                    curr_avail_paths = intersect(prev_avail_freqs_idx, find(full_costs(:, j) > 0));

                    min_cost_path = curr_avail_paths(1);
                    if (length(curr_avail_paths) >= 2) then
                        for (k = curr_avail_paths(2 : $))
                            if (cost0(k) < cost0(min_cost_path)) then
                                min_cost_path = k;
                            end
                        end
                    end
                    curr_total_full_cost(j) = cost0(min_cost_path);
                    opt_freqs_idx(frame_idx, j) = min_cost_path;
                end

                normalize_factor = max(curr_total_full_cost);
                curr_total_full_cost = curr_total_full_cost ./ normalize_factor;

                prev_avail_freqs_idx = curr_avail_freqs_idx;
                prev_total_full_cost = curr_total_full_cost;
                prev_freq_cndds = curr_freq_cndds;
                last_avail_freqs_idx = curr_avail_freqs_idx;
                last_valid_frm_idx = frame_idx;
                last_avl_tot_full_cost = curr_total_full_cost;
                is_prev_frm_valid = 1;
            else
                error('Internal fault.');
                return;
        end

        frame_idx = frame_idx + 1;
    end

    num_of_frames = frame_idx - 1;

    term_opt_freq_idx = 1;
    if (~isempty(last_avail_freqs_idx)) then
        term_opt_freq_idx = last_avail_freqs_idx(1);
        if (length(last_avail_freqs_idx) >= 2) then
            for k = last_avail_freqs_idx(2 : $)
                if (last_avl_tot_full_cost(k) < last_avl_tot_full_cost(term_opt_freq_idx)) then
                    term_opt_freq_idx = k;
                end
            end
        end
    else
        error('Base pitch detection failed. Increase frame_size or/and enlarge frequency search range [fmin, fmax], and then try again.');
        return;
    end

    y = [freq_candidates(last_valid_frm_idx, term_opt_freq_idx), zeros(1, num_of_frames - last_valid_frm_idx)];
    prev_opt_freq_idx = term_opt_freq_idx;
    q = 0;
    if (~isempty(part_term_frms_idx)) then
        if (last_valid_frm_idx < num_of_frames) then
            q = 2;
        else
            q = 1;
        end
    end

    for k = (last_valid_frm_idx - 1) : -1 : 1
        curr_opt_freq_idx = opt_freqs_idx(k + 1, prev_opt_freq_idx);
        if (curr_opt_freq_idx > 0) then
            y = [freq_candidates(k, curr_opt_freq_idx), y];
            prev_opt_freq_idx = curr_opt_freq_idx;
        else
            if ((q >= 1) & (q <= length(part_term_frms_idx))) then
                if ((k+1 == part_term_frms_idx(q)) & (k >= 2)) then
                    y = [freq_candidates(k - 1, part_opt_freqs_idx(q)), y];
                    prev_opt_freq_idx = part_opt_freqs_idx(q);
                    q = q + 1;
                else
                    y = [0, y];
                end
            else
                y = [0, y];
            end
        end
    end

    sz = size(y);
    t = linspace((n_start - 1) / fs, (n_end - 1) / fs, sz(2));
endfunction


