function [y, yf] = autbx_lpfilter(x, order, window_type, fs, cutoff_freq, fpar, y0)
    fc = cutoff_freq / fs;
    wft = wfir('lp', order, [fc, 0], window_type, fpar);

    rhs = argn(2);
    select (rhs)
    case 7 then
        [y, yf] = filter(wft, 1, x, y0);
    case 6 then
        [y, yf] = filter(wft, 1, x);
    else
        error('Wrong number of input arguments.');
        return;
    end
endfunction
