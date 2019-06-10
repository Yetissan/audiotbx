function [y] = autbx_medfilter(x, frame_size)
    y = [];

    x = (x(:))';
    N = length(x);

    // Mirrored extension on both ends.
    s = [x(frame_size/2:-1:1), x, x($:-1:$-frame_size/2+1)];

    for i = 1 : N
        y(i) = median(s(i:i+frame_size-1));
    end

    y = (y(:))';
endfunction
