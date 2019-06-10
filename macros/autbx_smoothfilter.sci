function [y, avg] = autbx_smoothfilter(x, num_of_classes, bounding_mode, du, dl, medfilt_order, out_mode)
    if ((num_of_classes < 2)) then
        error('num_of_classes must be greater than 2.');
        return;
    end

    x = (x(:))';
    if (medfilt_order >= 2) then
        x = autbx_medfilter(x, medfilt_order);
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
    
    select (bounding_mode)
    case 1 then
        ubound = avg + du;
        lbound = avg - dl;
    case 2 then
        ubound = max(v) * sqrt(2);
        lbound = min(v) ./ sqrt(2);
    else
        error('Invalid bounding mode. ');
        return;
    end

    u = [];
    for j = 1:length(x)
        if ((x(j) < ubound) & (x(j) > lbound)) then
            u = [u, [j; x(j)]];
        end
    end

    d = splin(u(1, :), u(2, :));
    y = interp(1:length(x), u(1, :), u(2, :), d, out_mode);
endfunction

