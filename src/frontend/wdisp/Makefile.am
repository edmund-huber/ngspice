## Process this file with automake to produce Makefile.in

noinst_LTLIBRARIES = libwindisp.la

libwindisp_la_SOURCES = \
	windisp.c	\
	windisp.h	\
	winprint.c	\
	winprint.h	\
	ftegraf.h


AM_CPPFLAGS = @AM_CPPFLAGS@ -I$(top_srcdir)/src/include -I$(top_srcdir)/src/frontend/wdisp

MAINTAINERCLEANFILES = Makefile.in
