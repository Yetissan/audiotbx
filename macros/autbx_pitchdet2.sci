function [y, t, avg] = autbx_pitchdet2(x, n_start, n_end, fs, n0, t1, t2, ...
        num_of_classes, bounding_mode, du, dl, medfilt_order)
    y = [];
    t = [];

    x = x(n_start:n_end);
    x = (x(:))';

    // First order peaks of selected portion of x.
    [xu, xl] = envelope(1:length(x), x, 1, 1);
    
    // Second order peaks of selected portion of x.
    [yu, yl] = envelope(xl(1, :), xl(2, :), 1, 1);
    xx = xl;
    yy = yl;
    
    if (length(yy(1, :)) == 0) then
        error('Unknown pitch.');
        return;
    end

    // Let yy(:, 1) be the first pitch marking point.
    n = yy(1, 1);
    u = yy(2, 1);
    j = 1;
    prev_len = fs / 2;

    for i = 2:length(yy(1, :))
        // Time interval between current pitch marking candidate and latest pitch marking
        // point.
        len = (yy(1, i) - n(j)) / fs;
        // First order peaks between current pitch marking candidate and latest pitch
        // marking point.
        fop = xx(:, find((xx(1, :) > n(j)) & (xx(1, :) < yy(1, i))));
        fopcnt = length(fop(1, :));

        p = find(fop(2, :) < yy(2, i));

        if ((len > t1) & (fopcnt >= n0) & isempty(p)) then
            n = [n, yy(1, i)];
            u = [u, yy(2, i)];
            j = j + 1;
            prev_len = len;
        elseif ((abs(len - prev_len) < prev_len/10) & (abs(yy(2,i) - u(j)) < u(j)/2)) then
            n = [n, yy(1, i)];
            u = [u, yy(2, i)];
            j = j + 1;
            prev_len = len;
        elseif ((len > t2) & (fopcnt > 0)) then
            for k = 1 : fopcnt
                len2 = (fop(1, k) - n(j)) / fs;
                if (abs(len2 - prev_len) < prev_len / 10) then
                    n = [n, fop(1, k)];
                    u = [u, fop(2, k)];
                    j = j + 1;
                    prev_len = len2;
                end
            end
        end
    end

    p = abs(n - [n(2:$), 0]);
    p($) = p($ - 1);
    [p, avg] = autbx_smoothfilter(p, num_of_classes, bounding_mode, du, dl, medfilt_order, 'natural');

    y = fs ./ p;
    z = size(y);
    t = linspace((n_start - 1)/fs, (n_end - 1)/fs, z(2));
endfunction
