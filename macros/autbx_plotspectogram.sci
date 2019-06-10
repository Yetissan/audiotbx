function autbx_plotspectogram(y, t, f, granularity)
    scf();
    grayplot(t, f, y);

    f = gcf();
    f.color_map = jetcolormap(granularity);

    a = gca();
    a.x_label.text = 'Time (s)';
    a.y_label.text = 'Frequency (Hz)';
    xtitle('Spectogram');
endfunction