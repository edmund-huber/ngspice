## Process this file with automake to produce Makefile.in

noinst_LTLIBRARIES = libtxl.la

libtxl_la_SOURCES = \
	txl.c      \
	txlacct.c  \
	txlask.c   \
	txldefs.h  \
	txldest.c  \
	txlext.h   \
	txlfbr.c   \
	txlinit.h  \
	txlitf.h   \
	txlload.c  \
	txlmdel.c  \
	txlparam.c \
	txldel.c   \
	txlmask.c  \
	txlmpar.c  \
	txlsetup.c \
	txlinit.c

AM_CPPFLAGS = @AM_CPPFLAGS@ -I$(top_srcdir)/src/include

MAINTAINERCLEANFILES = Makefile.in
