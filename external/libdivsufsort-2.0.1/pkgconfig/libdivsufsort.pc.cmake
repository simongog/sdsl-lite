prefix=@prefix@
exec_prefix=${prefix}
libdir=${exec_prefix}/lib
includedir=${prefix}/include

Name: @PROJECT_NAME@@W64BIT@
Description: @PROJECT_DESCRIPTION@
Version: @PROJECT_VERSION_FULL@
URL: @PROJECT_URL@
Libs: -L${libdir} -ldivsufsort@W64BIT@
Cflags: -I${includedir}
