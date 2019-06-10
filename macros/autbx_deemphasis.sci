function [y, yf] = autbx_deemphasis(x, mu, y0)
    rhs = argn(2);
    
    select (rhs)
    case 3 then
        [y, yf] = filter(1, [1, -mu], x, y0);
    case 2 then
        [y, yf] = filter(1, [1, -mu], x);
    else
        error('Wrong number of input arguments.');
    end
endfunction
