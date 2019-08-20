function builder_gw_c()

    includes_src_c = ilib_include_flag(get_absolute_file_path("builder_gateway_c.sce") + "../../src/c");

    WITHOUT_AUTO_PUTLHSVAR = %t;

    gateway_c_path = get_absolute_file_path("builder_gateway_c.sce");
    fcnsrctab = [   ...
//        "c_skeleton",       "sci_skeleton",             "sci_gwskeleton.c",
        "autbx_emd",        "sci_autbx_emd",            "sci_gw_autbx_emd.c"      ];

    tbx_build_gateway("audiotbx_c", fcnsrctab(:,1:2), fcnsrctab(:,3), gateway_c_path, ...
        ["../../src/c/" + LIB_NAME], "", includes_src_c);

endfunction

builder_gw_c();
clear builder_gw_c; // remove builder_gw_c on stack
