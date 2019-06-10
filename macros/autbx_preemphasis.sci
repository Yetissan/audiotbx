function [y, yf] = autbx_preemphasis(x, mu, y0)
    rhs = argn(2);
    
    select (rhs)
    case 3 then
        [y, yf] = filter([1, -mu], 1, x, y0);
    case 2 then
        [y, yf] = filter([1, -mu], 1, x);
    else
        error('Wrong number of input arguments.');
    end
endfunction
