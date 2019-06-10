function [y, yf] = autbx_bpfilter(x, order, window_type, fs, hp_cutoff_freq, lp_cutoff_freq, fpar, y0)
    fc = [hp_cutoff_freq, lp_cutoff_freq] / fs;
    wft = wfir('bp', order, fc, window_type, fpar);

    rhs = argn(2);
    select (rhs)
    case 8 then
        [y, yf] = filter(wft, 1, x, y0);
    case 7 then
        [y, yf] = filter(wft, 1, x);
    else
        error('Wrong number of input arguments.');
        return;
    end
endfunction
