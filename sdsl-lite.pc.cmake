prefix=@CMAKE_INSTALL_PREFIX@
exec_prefix=${prefix}
libdir=@CMAKE_INSTALL_LIBDIR@
includedir=@CMAKE_INSTALL_INCLUDEDIR@

Name: @PROJECT_NAME@
Description: @PROJECT_DESCRIPTION@
Version: @PROJECT_VERSION_FULL@
URL: @PROJECT_URL@
Libs: -L${libdir} -lsdsl -ldivsufsort -ldivsufsort64
Cflags: -I${includedir}
