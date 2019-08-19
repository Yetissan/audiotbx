// This macro compiles the files

function builder_c()
    src_c_path = get_absolute_file_path("builder_c.sce");

    CFLAGS = ilib_include_flag(src_c_path);

    fcnsrctab = [ ...
        "EmpiricalModeDecomposition", "EmpiricalModeDecomposition.c" ];

    tbx_build_src(fcnsrctab(:,1), fcnsrctab(:,2), "c", src_c_path, "", ...
        "", CFLAGS, "", "", LIB_NAME);
endfunction

builder_c();
clear builder_c; // remove builder_c on stack

