# Makefile.am for libdivsufsort

lib_LTLIBRARIES = libdivsufsort.la
libdivsufsort_la_SOURCES = divsufsort.c sssort.c trsort.c utils.c
libdivsufsort_la_CPPFLAGS = -I$(top_srcdir)/include
libdivsufsort_la_LDFLAGS = -version-info @LT_CURRENT@:@LT_REVISION@:@LT_AGE@ -export-symbols $(srcdir)/libdivsufsort.sym

if DIVSUFSORT64
lib_LTLIBRARIES += libdivsufsort64.la
libdivsufsort64_la_SOURCES = divsufsort.c sssort.c trsort.c utils.c
libdivsufsort64_la_CPPFLAGS = -I$(top_srcdir)/include -DBUILD_DIVSUFSORT64
libdivsufsort64_la_LDFLAGS = -version-info @LT_CURRENT@:@LT_REVISION@:@LT_AGE@ -export-symbols $(srcdir)/libdivsufsort64.sym
endif

EXTRA_DIST = libdivsufsort.sym libdivsufsort64.sym CMakeLists.txt
