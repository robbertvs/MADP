# $Id: Makefile.am 9370 2015-03-10 10:42:38Z mtjspaan $

# Lists all the subdirs where 'make' has to look
# When changing the subdirs, also change AC_CONFIG_FILES in configure.ac
SUBDIRS_LIBS = libDAI base support parser planning
if !DARWIN
SUBDIRS_LIBS += libpomdp-solve
endif
SUBDIRS_BINARIES = utils examples solvers 
# first the component libraries, then the madp to build the shared
# library followed by the binaries
SUBDIRS = $(SUBDIRS_LIBS) madp $(SUBDIRS_BINARIES)

include $(top_srcdir)/src/Makefile.custom
include $(top_srcdir)/config/Make-rules

generaldir = $(pkgincludedir)
general_HEADERS = include/versions.h include/configuration.h

htmldoc:
	cd ../doc; doxygen doxygen.conf

