## Process this file with automake to produce Makefile.in

noinst_LTLIBRARIES = libmathmisc.la

libmathmisc_la_SOURCES = \
	accuracy.c	\
	accuracy.h	\
	bernoull.h	\
	bernoull.c	\
	erfc.c		\
	equality.c	\
	isinf.c		\
	isnan.c		\
	logb.c		\
	scalb.c		\
	norm.h		\
	norm.c		\
	randnumb.c


EXTRA_DIST = test_accuracy.c test_erfc.c

AM_CPPFLAGS = @AM_CPPFLAGS@ -I$(top_srcdir)/src/include

MAINTAINERCLEANFILES = Makefile.in
