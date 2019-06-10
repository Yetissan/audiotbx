function builder_src()
    language_src = ["c"];
    path_src = get_absolute_file_path("builder_src.sce");
    tbx_builder_src_lang(language_src, path_src);
endfunction

builder_src();
clear builder_src; // remove builder_src on stack
