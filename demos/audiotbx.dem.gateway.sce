function subdemolist = demo_gateway()
    demopath = get_absolute_file_path("audiotbx.dem.gateway.sce");

// -+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
    subdemolist = [ ...
        "Skeleton of a demo",                   "skeleton.dem.sce"];

    subdemolist(:,2) = demopath + subdemolist(:,2);

endfunction

subdemolist = demo_gateway();
clear demo_gateway; // remove demo_gateway on stack
