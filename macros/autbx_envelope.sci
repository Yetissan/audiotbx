function [yu, yl] = autbx_envelope(pos, x, strict)
    yu = [];
    yl = [];
    N = length(x);
    
    if (N < 4) then
        error('Not enough data provided.');
        return;
    end
    
    if (length(pos) ~= N) then
        error('Vector pos and x must have equal length.');
        return;
    end
    
    if (x(1) >= x(2)) then
        yu = [yu, [pos(1); x(1)]];
    end

    if (x(1) <= x(2)) then
        yl = [yl, [pos(1); x(1)]];
    end
    
    if (strict == 1) then
        for i = 2:N-1
            if (((x(i) >= x(i-1)) & (x(i) > x(i+1))) | ((x(i) > x(i-1)) & (x(i) >= x(i+1)))) then
                yu = [yu, [pos(i); x(i)]];
            end
    
            if (((x(i) <= x(i-1)) & (x(i) < x(i+1))) | ((x(i) < x(i-1)) & (x(i) <= x(i+1)))) then
                yl = [yl, [pos(i); x(i)]];
            end
        end
    else
        for i = 2:N-1
            if ((x(i) >= x(i-1)) & (x(i) >= x(i+1))) then
                yu = [yu, [pos(i); x(i)]];
            end
    
            if ((x(i) <= x(i-1)) & (x(i) <= x(i+1))) then
                yl = [yl, [pos(i); x(i)]];
            end
        end
    end

    if (x(N) >= x(N-1)) then
        yu = [yu, [pos(N); x(N)]];
    end

    if (x(N) <= x(N-1)) then
        yl = [yl, [pos(N); x(N)]];
    end
endfunction