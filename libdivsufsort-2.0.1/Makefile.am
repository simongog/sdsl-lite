# Makefile.am for libdivsufsort

SUBDIRS = include lib examples

EXTRA_DIST = ChangeLog.old CMakeLists.txt VERSION CMakeModules pkgconfig

ACLOCAL_AMFLAGS = -I m4

libtool: $(LIBTOOL_DEPS)
	$(SHELL) ./config.status --recheck
