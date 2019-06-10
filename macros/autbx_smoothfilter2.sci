function [y, avg] = smooth_filter2(x, num_of_classes, th1, th2)
    x = (x(:))';
    
    if (length(x) <= 5) then
        y = x;
        return;
    end
    
    if ((num_of_classes < 2)) then
        error('num_of_classes must be greater than 2.');
        return;
    end

    n = linspace(min(x), max(x), num_of_classes);
    // printf('min(x) = %g, max(x) = %g, num_of_classes = %d.\n', min(x), max(x), num_of_classes);
    histdata = [];
    for i = 2:length(n)
        histdata = [histdata, length(find((x > n(i-1)) & (x <= n(i))))];
    end

    N = length(histdata);
    m = find(histdata >= (0.8.*max(histdata)));
    // Find first apperance of contigeous subsequence from m, assuming that x is unimodal
    // distributed.
    u = m(1);
    if (length(m) >= 2) then
        for i = 2 : length(m)
            if (m(i) == u(i-1)+1) then
                u = [u, m(i)];
            else
                break;
            end
        end
    end

    v = ((max(x) - min(x))./ (num_of_classes - 1)) .* u + min(x);
    avg = sum(v) ./ length(v);
    
    ubound = avg .* sqrt(2);
    lbound = avg ./ sqrt(2);
    u = [];
    for j = 1:length(x)
        if ((x(j) < ubound) & (x(j) > lbound)) then
            u = [u, [j; x(j)]];
        end
    end

    d = splin(u(1, :), u(2, :));
    y = interp(1:length(x)+2, u(1, :), u(2, :), d, 'natural');
    
    for j = 2 : (length(y) - 1)
        if (abs(y(j) ./ 2 - y(j-1)) < th1) then
            y(j) = y(j) ./ 2; 
        elseif (abs(2 .* y(j) - y(j-1)) < th1) then
            y(j) = 2 .* y(j);
        end
        
        if ((abs(y(j) - y(j-1)) > th1) & (abs(y(j+1) - y(j-1)) > th2)) then
            if (j > 2) then
                y(j) = 2 .* y(j-1) - y(j-2);
            else
                y(j) = y(j-1);
            end
        elseif ((abs(y(j) - y(j-1)) > th1) & (abs(y(j+1) - y(j-1)) < th2)) then
            y(j) = (y(j-1) + y(j+1)) ./ 2;
        end
    end

    y = y(2:$-1);
endfunction
