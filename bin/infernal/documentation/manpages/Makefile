top_srcdir = ../..
srcdir     = .

SHELL      = /bin/sh

# Installation targets
#
prefix      = /usr/local
exec_prefix = ${prefix}
datarootdir = ${prefix}/share
bindir      = ${exec_prefix}/bin
libdir      = ${exec_prefix}/lib
includedir  = ${prefix}/include
mandir      = ${datarootdir}/man
docdir      = ${datarootdir}/doc/${PACKAGE_TARNAME}
pdfdir      = ${docdir}
mandir      = ${datarootdir}/man
man1dir     = ${mandir}/man1
man1ext     = .1

INSTALL     = /usr/bin/install -c

MANS =  infernal\
	cmalign\
	cmbuild\
	cmcalibrate\
	cmconvert\
	cmemit\
	cmfetch\
	cmpress\
	cmsearch\
	cmscan\
	cmstat
# don't install these, since the programs aren't installed
#	sse_cmsearch\
#	cmsearch_dual

.PHONY:  install uninstall clean distclean

install:
	for file in ${MANS}; do \
	   ${INSTALL} -m 0755 ${srcdir}/$$file.man ${DESTDIR}${man1dir}/$${file}${man1ext} ;\
	done

uninstall:
	for file in ${MANS}; do \
	   rm ${DESTDIR}${man1dir}/$${file}${man1ext} ;\
	done

clean:
	-rm -f *~ 

distclean:
	-rm -f *~ 
	-rm -f Makefile

################################################################
# Infernal - inference of RNA secondary structure alignments
# Version 1.1.1; July 2014
# Copyright (C) 2014 Howard Hughes Medical Institute.
# Other copyrights also apply. See the COPYRIGHT file for a full list.
# 
# Infernal is distributed under the terms of the GNU General Public License
# (GPLv3). See the LICENSE file for details.
################################################################
