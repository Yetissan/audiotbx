funcprot(0);

function [y] = autbx_viterbi(s, transcost, initcost)
    [num_of_obsvs, num_of_states] = size(s);

    if (num_of_obsvs <= 1) then
        error('Number of obervations must be greater or equal to 2.');
        return;
    end

    if (num_of_states <= 1) then
        error('Number of states must be greater or equal to 2.');
        return;
    end

    // Initialization
    opt_nodes = zeros(num_of_obsvs, num_of_states);
    prev_cost = initcost;
    prev_cost = (prev_cost(:))';
    prev_avail_nodes = find(s(1, :) > 0);
    prev_avail_nodes = (prev_avail_nodes(:))';
    if (isempty(prev_avail_nodes)) then
        error('An empty observation vector was found.');
        return;
    end
    
    // Recursion
    for i = 2 : num_of_obsvs
        A = transcost(i);
        s0 = s(i, :);
        curr_avail_nodes = find(s0 > 0);
        curr_avail_nodes = (curr_avail_nodes(:))';
        if (isempty(curr_avail_nodes)) then
            error('An empty observation vector was found.');
            return;
        end
        
        cost = zeros(1, num_of_states);
        for j = curr_avail_nodes
            cost0 = prev_cost .* A(:, j)';
            avail_paths = intersect(find(A(:, j) > 0), prev_avail_nodes);
            if (isempty(avail_paths)) then
                error('Broken path.');
                return;
            end
            min_cost_path = avail_paths(1);
            if (length(avail_paths) >= 2) then
                for k = avail_paths(2:$)
                    if (cost0(k) < cost0(min_cost_path)) then
                        min_cost_path = k;
                    end
                end
            end
            cost(j) = cost0(min_cost_path);
            opt_nodes(i, j) = min_cost_path;
        end

        prev_avail_nodes = curr_avail_nodes;
        prev_cost = cost(:)';
    end

    // Termination
    term_opt_node = prev_avail_nodes(1);
    if (length(prev_avail_nodes) >= 2) then
        for k = prev_avail_nodes(2:$)
            if (prev_cost(k) < prev_cost(term_opt_node)) then
                term_opt_node = k;
            end
        end
    end

    // Backtracking
    y = term_opt_node;
    prev_opt_node = term_opt_node;
    for l = num_of_obsvs - 1 : -1 : 1
        curr_opt_node = opt_nodes(l+1, prev_opt_node);
        y = [curr_opt_node, y];
        prev_opt_node = curr_opt_node;
    end
endfunction


