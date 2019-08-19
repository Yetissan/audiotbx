// This file is released under the 3-clause BSD license. See COPYING-BSD.
// Generated by builder.sce : Please, do not edit this file
// ----------------------------------------------------------------------------
//
if ~win64() then
  warning(_("This module requires a Windows x64 platform."));
  return
end
//
audiotbx_c_path = get_absolute_file_path('loader.sce');
//
// ulink previous function with same name
[bOK, ilib] = c_link('audiotbx_c');
if bOK then
  ulink(ilib);
end
//
link(audiotbx_c_path + filesep() + '../../src/c/libaudiotbx' + getdynlibext());
list_functions = [ 'autbx_emd';
];
addinter(audiotbx_c_path + filesep() + 'audiotbx_c' + getdynlibext(), 'audiotbx_c', list_functions);
// remove temp. variables on stack
clear audiotbx_c_path;
clear bOK;
clear ilib;
clear list_functions;
// ----------------------------------------------------------------------------
