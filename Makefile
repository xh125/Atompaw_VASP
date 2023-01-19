# Makefile.in generated by automake 1.12.4 from Makefile.am.
# src/Makefile.  Generated from Makefile.in by configure.

# Copyright (C) 1994-2012 Free Software Foundation, Inc.

# This Makefile.in is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, to the extent permitted by law; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.



# -*- Automake -*-
#
# Makefile for the AtomPAW package
#
# Copyright (C) 2010 Yann Pouillon
#
# This file is part of the AtomPAW software package. For license information,
# please see the COPYING file in the top-level directory of the source
# distribution.
#



am__make_dryrun = \
  { \
    am__dry=no; \
    case $$MAKEFLAGS in \
      *\\[\ \	]*) \
        echo 'am--echo: ; @echo "AM"  OK' | $(MAKE) -f - 2>/dev/null \
          | grep '^AM OK$$' >/dev/null || am__dry=yes;; \
      *) \
        for am__flg in $$MAKEFLAGS; do \
          case $$am__flg in \
            *=*|--*) ;; \
            *n*) am__dry=yes; break;; \
          esac; \
        done;; \
    esac; \
    test $$am__dry = yes; \
  }
pkgdatadir = $(datadir)/atompaw
pkgincludedir = $(includedir)/atompaw
pkglibdir = $(libdir)/atompaw
pkglibexecdir = $(libexecdir)/atompaw
am__cd = CDPATH="$${ZSH_VERSION+.}$(PATH_SEPARATOR)" && cd
install_sh_DATA = $(install_sh) -c -m 644
install_sh_PROGRAM = $(install_sh) -c
install_sh_SCRIPT = $(install_sh) -c
INSTALL_HEADER = $(INSTALL_DATA)
transform = $(program_transform_name)
NORMAL_INSTALL = :
PRE_INSTALL = :
POST_INSTALL = :
NORMAL_UNINSTALL = :
PRE_UNINSTALL = :
POST_UNINSTALL = :
build_triplet = x86_64-unknown-linux-gnu
host_triplet = x86_64-unknown-linux-gnu
target_triplet = x86_64-unknown-linux-gnu
bin_PROGRAMS = atompaw$(EXEEXT) graphatom$(EXEEXT)
subdir = src
DIST_COMMON = $(srcdir)/Makefile.am $(srcdir)/Makefile.in \
	$(srcdir)/pkginfo.f90.in
ACLOCAL_M4 = $(top_srcdir)/aclocal.m4
am__aclocal_m4_deps =  \
	$(top_srcdir)/config/m4/ax_f90_module_extension.m4 \
	$(top_srcdir)/config/m4/fortran.m4 \
	$(top_srcdir)/config/m4/libtool.m4 \
	$(top_srcdir)/config/m4/libxc.m4 \
	$(top_srcdir)/config/m4/libxc_names.m4 \
	$(top_srcdir)/config/m4/linalg.m4 \
	$(top_srcdir)/config/m4/ltoptions.m4 \
	$(top_srcdir)/config/m4/ltsugar.m4 \
	$(top_srcdir)/config/m4/ltversion.m4 \
	$(top_srcdir)/config/m4/lt~obsolete.m4 \
	$(top_srcdir)/config/m4/optimizations.m4 \
	$(top_srcdir)/config/m4/workarounds.m4 \
	$(top_srcdir)/configure.ac
am__configure_deps = $(am__aclocal_m4_deps) $(CONFIGURE_DEPENDENCIES) \
	$(ACLOCAL_M4)
mkinstalldirs = $(install_sh) -d
CONFIG_HEADER = $(top_builddir)/config.h
CONFIG_CLEAN_FILES = pkginfo.f90
CONFIG_CLEAN_VPATH_FILES =
am__vpath_adj_setup = srcdirstrip=`echo "$(srcdir)" | sed 's|.|.|g'`;
am__vpath_adj = case $$p in \
    $(srcdir)/*) f=`echo "$$p" | sed "s|^$$srcdirstrip/||"`;; \
    *) f=$$p;; \
  esac;
am__strip_dir = f=`echo $$p | sed -e 's|^.*/||'`;
am__install_max = 40
am__nobase_strip_setup = \
  srcdirstrip=`echo "$(srcdir)" | sed 's/[].[^$$\\*|]/\\\\&/g'`
am__nobase_strip = \
  for p in $$list; do echo "$$p"; done | sed -e "s|$$srcdirstrip/||"
am__nobase_list = $(am__nobase_strip_setup); \
  for p in $$list; do echo "$$p $$p"; done | \
  sed "s| $$srcdirstrip/| |;"' / .*\//!s/ .*/ ./; s,\( .*\)/[^/]*$$,\1,' | \
  $(AWK) 'BEGIN { files["."] = "" } { files[$$2] = files[$$2] " " $$1; \
    if (++n[$$2] == $(am__install_max)) \
      { print $$2, files[$$2]; n[$$2] = 0; files[$$2] = "" } } \
    END { for (dir in files) print dir, files[dir] }'
am__base_list = \
  sed '$$!N;$$!N;$$!N;$$!N;$$!N;$$!N;$$!N;s/\n/ /g' | \
  sed '$$!N;$$!N;$$!N;$$!N;s/\n/ /g'
am__uninstall_files_from_dir = { \
  test -z "$$files" \
    || { test ! -d "$$dir" && test ! -f "$$dir" && test ! -r "$$dir"; } \
    || { echo " ( cd '$$dir' && rm -f" $$files ")"; \
         $(am__cd) "$$dir" && rm -f $$files; }; \
  }
am__installdirs = "$(DESTDIR)$(libdir)" "$(DESTDIR)$(bindir)"
LTLIBRARIES = $(lib_LTLIBRARIES)
libatompaw_la_LIBADD =
am__objects_1 = abinitinterface.lo aeatom.lo anderson_driver.lo \
	atomdata.lo atompaw_report.lo blockdavidson_mod.lo excor.lo \
	exx_mod.lo exx_pseudo.lo exxdata.lo fock.lo general_mod.lo \
	globalmath.lo graphatom_report.lo gridmod.lo hf_mod.lo \
	hf_pseudo.lo interpolation_mod.lo ldagga_mod.lo libxc_mod.lo \
	numerov_mod.lo paw_sub.lo pseudo_sub.lo pseudo.lo \
	pseudodata.lo pwscfinterface.lo radialsr.lo report_mod.lo \
	search_sort.lo xmlinterface.lo \
	ini.lo asa.lo optreal.lo radial.lo setlocalpp.lo \
	xclib.lo xclib_grad.lo \
	pseudo_vasp.lo  cl_shift.lo rhfatm.lo  vasp_pseudo.lo 
am_libatompaw_la_OBJECTS = $(am__objects_1) atompaw_lib.lo
am__objects_2 = pkginfo.lo
nodist_libatompaw_la_OBJECTS = $(am__objects_2)
libatompaw_la_OBJECTS = $(am_libatompaw_la_OBJECTS) \
	$(nodist_libatompaw_la_OBJECTS)
libatompaw_la_LINK = $(LIBTOOL) --tag=FC $(AM_LIBTOOLFLAGS) \
	$(LIBTOOLFLAGS) --mode=link $(FCLD) $(AM_FCFLAGS) $(FCFLAGS) \
	$(libatompaw_la_LDFLAGS) $(LDFLAGS) -o $@
PROGRAMS = $(bin_PROGRAMS)
am_atompaw_OBJECTS = atompaw_prog.$(OBJEXT)
atompaw_OBJECTS = $(am_atompaw_OBJECTS)
atompaw_DEPENDENCIES = libatompaw.la
am_graphatom_OBJECTS = graphatom.$(OBJEXT)
graphatom_OBJECTS = $(am_graphatom_OBJECTS)
graphatom_DEPENDENCIES = libatompaw.la
DEFAULT_INCLUDES = -I. -I$(top_builddir)
PPFCCOMPILE = $(FC) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
	$(AM_CPPFLAGS) $(CPPFLAGS) $(AM_FCFLAGS) $(FCFLAGS)
LTPPFCCOMPILE = $(LIBTOOL) --tag=FC $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) \
	--mode=compile $(FC) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
	$(AM_CPPFLAGS) $(CPPFLAGS) $(AM_FCFLAGS) $(FCFLAGS)
FCLD = $(FC)
FCLINK = $(LIBTOOL) --tag=FC $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) \
	--mode=link $(FCLD) $(AM_FCFLAGS) $(FCFLAGS) $(AM_LDFLAGS) \
	$(LDFLAGS) -o $@
FCCOMPILE = $(FC) $(AM_FCFLAGS) $(FCFLAGS)
LTFCCOMPILE = $(LIBTOOL) --tag=FC $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) \
	--mode=compile $(FC) $(AM_FCFLAGS) $(FCFLAGS)
SOURCES = $(libatompaw_la_SOURCES) $(nodist_libatompaw_la_SOURCES) \
	$(atompaw_SOURCES) $(graphatom_SOURCES)
DIST_SOURCES = $(libatompaw_la_SOURCES) $(atompaw_SOURCES) \
	$(graphatom_SOURCES)
am__can_run_installinfo = \
  case $$AM_UPDATE_INFO_DIR in \
    n|no|NO) false;; \
    *) (install-info --version) >/dev/null 2>&1;; \
  esac
ETAGS = etags
CTAGS = ctags
DISTFILES = $(DIST_COMMON) $(DIST_SOURCES) $(TEXINFOS) $(EXTRA_DIST)
ACLOCAL = ${SHELL} /home/jun_jiang/Softs/atompaw_test/config/gnu/missing --run aclocal-1.12
AMTAR = $${TAR-tar}
AR = ar
ATOMPAW_LIBDIR = 
ATP_FCOPTS = -O2
ATP_LDOPTS = -O2 -Vaxlib
ATP_LIBS = 
AUTOCONF = ${SHELL} /home/jun_jiang/Softs/atompaw_test/config/gnu/missing --run autoconf
AUTOHEADER = ${SHELL} /home/jun_jiang/Softs/atompaw_test/config/gnu/missing --run autoheader
AUTOMAKE = ${SHELL} /home/jun_jiang/Softs/atompaw_test/config/gnu/missing --run automake-1.12
AWK = gawk
CC = icc
CCDEPMODE = depmode=none
CFLAGS = -g -O2
CPP = icc -E
CPPFLAGS = 
CYGPATH_W = echo
DEFS = -DHAVE_CONFIG_H
DEPDIR = .deps
DLLTOOL = false
DSYMUTIL = 
DUMPBIN = 
DVIPDF = dvipdf
ECHO_C = 
ECHO_N = -n
ECHO_T = 
EGREP = /bin/grep -E
EXEEXT = 
FC = ifort
FCFLAGS = -g
FCFLAGS_F90 = 
FGREP = /bin/grep -F
GREP = /bin/grep
INSTALL = /usr/bin/install -c
INSTALL_DATA = ${INSTALL} -m 644
INSTALL_PROGRAM = ${INSTALL}
INSTALL_SCRIPT = ${INSTALL}
INSTALL_STRIP_PROGRAM = $(install_sh) -c -s
LATEX = latex
LD = /usr/bin/ld -m elf_x86_64
LDFLAGS = -O2 -Vaxlib
LIBOBJS = 
LIBS = -L/opt/intel/composer_xe_2013.0.079/mkl/lib/intel64 -lmkl_core -lmkl_intel_lp64 -lmkl_intel_thread -liomp5  
LIBTOOL = $(SHELL) $(top_builddir)/libtool
LIPO = 
LN_S = ln -s
LTLIBOBJS = 
LTOBJEXT = lo
MAKEINFO = ${SHELL} /home/jun_jiang/Softs/atompaw_test/config/gnu/missing --run makeinfo
MANIFEST_TOOL = :
MKDIR_P = /bin/mkdir -p
MODEXT = mod
NM = /usr/bin/nm -B
NMEDIT = 
OBJDUMP = objdump
OBJEXT = o
OTOOL = 
OTOOL64 = 
PACKAGE = atompaw
PACKAGE_BUGREPORT = natalie@wfu.edu, marc.torrent@cea.fr
PACKAGE_NAME = AtomPAW
PACKAGE_STRING = AtomPAW 4.0.0.4
PACKAGE_TARNAME = atompaw
PACKAGE_URL = 
PACKAGE_VERSION = 4.0.0.4
PATH_SEPARATOR = :
PDFLATEX = pdflatex
PERL = perl
RANLIB = ranlib
SED = /bin/sed
SET_MAKE = 
SHELL = /bin/bash
STRIP = strip
VERSION = 4.0.0.4
abs_builddir = /home/jun_jiang/Softs/atompaw_test/src
abs_srcdir = /home/jun_jiang/Softs/atompaw_test/src
abs_top_builddir = /home/jun_jiang/Softs/atompaw_test
abs_top_srcdir = /home/jun_jiang/Softs/atompaw_test
ac_ct_AR = ar
ac_ct_CC = icc
ac_ct_DUMPBIN = 
ac_ct_FC = 
am__include = include
am__leading_dot = .
am__quote = 
am__tar = $${TAR-tar} chof - "$$tardir"
am__untar = $${TAR-tar} xf -
atp_fc_path = /opt/intel/composer_xe_2013.0.079/bin/intel64/ifort
atp_fc_vendor = intel
atp_fc_version = 0.0
atp_fc_wrap = no
bindir = ${exec_prefix}/bin
build = x86_64-unknown-linux-gnu
build_alias = 
build_cpu = x86_64
build_os = linux-gnu
build_vendor = unknown
builddir = .
datadir = ${datarootdir}
datarootdir = ${prefix}/share
docdir = ${datarootdir}/doc/${PACKAGE_TARNAME}
dvidir = ${docdir}
exec_prefix = ${prefix}
host = x86_64-unknown-linux-gnu
host_alias = 
host_cpu = x86_64
host_os = linux-gnu
host_vendor = unknown
htmldir = ${docdir}
includedir = ${prefix}/include
infodir = ${datarootdir}/info
install_sh = ${SHELL} /home/jun_jiang/Softs/atompaw_test/config/gnu/install-sh
libdir = ${exec_prefix}/lib
libexecdir = ${exec_prefix}/libexec
localedir = ${datarootdir}/locale
localstatedir = ${prefix}/var
mandir = ${datarootdir}/man
mkdir_p = $(MKDIR_P)
oldincludedir = /usr/include
pdfdir = ${docdir}
prefix = /home/jun_jiang/Softs/atompaw_test
program_transform_name = s,x,x,
psdir = ${docdir}
sbindir = ${exec_prefix}/sbin
sharedstatedir = ${prefix}/com
srcdir = .
sysconfdir = ${prefix}/etc
target = x86_64-unknown-linux-gnu
target_alias = 
target_cpu = x86_64
target_os = linux-gnu
target_vendor = unknown
top_build_prefix = ../
top_builddir = ..
top_srcdir = ..

#
# Main source files
#

# Common source files
atp_srcs = \
  abinitinterface.f90 \
  aeatom.f90 \
  anderson_driver.f90 \
  atomdata.f90 \
  atompaw_report.f90 \
  blockdavidson_mod.f90 \
  excor.f90 \
  exx_mod.f90 \
  exx_pseudo.f90 \
  exxdata.f90 \
  fock.f90 \
  general_mod.f90 \
  globalmath.f90 \
  graphatom_report.f90 \
  gridmod.f90 \
  hf_mod.f90 \
  hf_pseudo.f90 \
  interpolation_mod.f90 \
  ldagga_mod.f90 \
  libxc_mod.F90 \
  numerov_mod.f90 \
  paw_sub.f90 \
  pseudo_sub.f90 \
  pseudo.f90 \
  ini.f90 \
  asa.f90 \
  optreal.f90 \
  radial.f90 \
  xclib.f90 \
  xclib_grad.f90 \
  setlocalpp.f90 \
  pseudo_vasp.f90 \
  cl_shift.f90 \
  rhfatm.f90 \
  vasp_pseudo.f90 \
  pseudodata.f90 \
  pwscfinterface.f90 \
  radialsr.f90 \
  report_mod.f90 \
  search_sort.f90 \
  xmlinterface.f90


# Built source files (generated by configure)
atp_built_srcs = \
  pkginfo.f90


# Libraries to install
lib_LTLIBRARIES = libatompaw.la
libatompaw_la_SOURCES = $(atp_srcs) atompaw_lib.f90
nodist_libatompaw_la_SOURCES = $(atp_built_srcs)
libatompaw_la_LDFLAGS = -version-info 0:0:0
atompaw_SOURCES = atompaw_prog.f90
atompaw_LDADD = libatompaw.la
graphatom_SOURCES = graphatom.f90
graphatom_LDADD = libatompaw.la
all: all-am

.SUFFIXES:
.SUFFIXES: .F90 .f90 .lo .o .obj
$(srcdir)/Makefile.in:  $(srcdir)/Makefile.am  $(am__configure_deps)
	@for dep in $?; do \
	  case '$(am__configure_deps)' in \
	    *$$dep*) \
	      ( cd $(top_builddir) && $(MAKE) $(AM_MAKEFLAGS) am--refresh ) \
	        && { if test -f $@; then exit 0; else break; fi; }; \
	      exit 1;; \
	  esac; \
	done; \
	echo ' cd $(top_srcdir) && $(AUTOMAKE) --gnu src/Makefile'; \
	$(am__cd) $(top_srcdir) && \
	  $(AUTOMAKE) --gnu src/Makefile
.PRECIOUS: Makefile
Makefile: $(srcdir)/Makefile.in $(top_builddir)/config.status
	@case '$?' in \
	  *config.status*) \
	    cd $(top_builddir) && $(MAKE) $(AM_MAKEFLAGS) am--refresh;; \
	  *) \
	    echo ' cd $(top_builddir) && $(SHELL) ./config.status $(subdir)/$@ $(am__depfiles_maybe)'; \
	    cd $(top_builddir) && $(SHELL) ./config.status $(subdir)/$@ $(am__depfiles_maybe);; \
	esac;

$(top_builddir)/config.status: $(top_srcdir)/configure $(CONFIG_STATUS_DEPENDENCIES)
	cd $(top_builddir) && $(MAKE) $(AM_MAKEFLAGS) am--refresh

$(top_srcdir)/configure:  $(am__configure_deps)
	cd $(top_builddir) && $(MAKE) $(AM_MAKEFLAGS) am--refresh
$(ACLOCAL_M4):  $(am__aclocal_m4_deps)
	cd $(top_builddir) && $(MAKE) $(AM_MAKEFLAGS) am--refresh
$(am__aclocal_m4_deps):
pkginfo.f90: $(top_builddir)/config.status $(srcdir)/pkginfo.f90.in
	cd $(top_builddir) && $(SHELL) ./config.status $(subdir)/$@
install-libLTLIBRARIES: $(lib_LTLIBRARIES)
	@$(NORMAL_INSTALL)
	@list='$(lib_LTLIBRARIES)'; test -n "$(libdir)" || list=; \
	list2=; for p in $$list; do \
	  if test -f $$p; then \
	    list2="$$list2 $$p"; \
	  else :; fi; \
	done; \
	test -z "$$list2" || { \
	  echo " $(MKDIR_P) '$(DESTDIR)$(libdir)'"; \
	  $(MKDIR_P) "$(DESTDIR)$(libdir)" || exit 1; \
	  echo " $(LIBTOOL) $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) --mode=install $(INSTALL) $(INSTALL_STRIP_FLAG) $$list2 '$(DESTDIR)$(libdir)'"; \
	  $(LIBTOOL) $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) --mode=install $(INSTALL) $(INSTALL_STRIP_FLAG) $$list2 "$(DESTDIR)$(libdir)"; \
	}

uninstall-libLTLIBRARIES:
	@$(NORMAL_UNINSTALL)
	@list='$(lib_LTLIBRARIES)'; test -n "$(libdir)" || list=; \
	for p in $$list; do \
	  $(am__strip_dir) \
	  echo " $(LIBTOOL) $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) --mode=uninstall rm -f '$(DESTDIR)$(libdir)/$$f'"; \
	  $(LIBTOOL) $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) --mode=uninstall rm -f "$(DESTDIR)$(libdir)/$$f"; \
	done

clean-libLTLIBRARIES:
	-test -z "$(lib_LTLIBRARIES)" || rm -f $(lib_LTLIBRARIES)
	@list='$(lib_LTLIBRARIES)'; \
	locs=`for p in $$list; do echo $$p; done | \
	      sed 's|^[^/]*$$|.|; s|/[^/]*$$||; s|$$|/so_locations|' | \
	      sort -u`; \
	test -z "$$locs" || { \
	  echo rm -f $${locs}; \
	  rm -f $${locs}; \
	}
libatompaw.la: $(libatompaw_la_OBJECTS) $(libatompaw_la_DEPENDENCIES) $(EXTRA_libatompaw_la_DEPENDENCIES) 
	$(libatompaw_la_LINK) -rpath $(libdir) $(libatompaw_la_OBJECTS) $(libatompaw_la_LIBADD) $(LIBS)
install-binPROGRAMS: $(bin_PROGRAMS)
	@$(NORMAL_INSTALL)
	@list='$(bin_PROGRAMS)'; test -n "$(bindir)" || list=; \
	if test -n "$$list"; then \
	  echo " $(MKDIR_P) '$(DESTDIR)$(bindir)'"; \
	  $(MKDIR_P) "$(DESTDIR)$(bindir)" || exit 1; \
	fi; \
	for p in $$list; do echo "$$p $$p"; done | \
	sed 's/$(EXEEXT)$$//' | \
	while read p p1; do if test -f $$p || test -f $$p1; \
	  then echo "$$p"; echo "$$p"; else :; fi; \
	done | \
	sed -e 'p;s,.*/,,;n;h' -e 's|.*|.|' \
	    -e 'p;x;s,.*/,,;s/$(EXEEXT)$$//;$(transform);s/$$/$(EXEEXT)/' | \
	sed 'N;N;N;s,\n, ,g' | \
	$(AWK) 'BEGIN { files["."] = ""; dirs["."] = 1 } \
	  { d=$$3; if (dirs[d] != 1) { print "d", d; dirs[d] = 1 } \
	    if ($$2 == $$4) files[d] = files[d] " " $$1; \
	    else { print "f", $$3 "/" $$4, $$1; } } \
	  END { for (d in files) print "f", d, files[d] }' | \
	while read type dir files; do \
	    if test "$$dir" = .; then dir=; else dir=/$$dir; fi; \
	    test -z "$$files" || { \
	    echo " $(INSTALL_PROGRAM_ENV) $(LIBTOOL) $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) --mode=install $(INSTALL_PROGRAM) $$files '$(DESTDIR)$(bindir)$$dir'"; \
	    $(INSTALL_PROGRAM_ENV) $(LIBTOOL) $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) --mode=install $(INSTALL_PROGRAM) $$files "$(DESTDIR)$(bindir)$$dir" || exit $$?; \
	    } \
	; done

uninstall-binPROGRAMS:
	@$(NORMAL_UNINSTALL)
	@list='$(bin_PROGRAMS)'; test -n "$(bindir)" || list=; \
	files=`for p in $$list; do echo "$$p"; done | \
	  sed -e 'h;s,^.*/,,;s/$(EXEEXT)$$//;$(transform)' \
	      -e 's/$$/$(EXEEXT)/' `; \
	test -n "$$list" || exit 0; \
	echo " ( cd '$(DESTDIR)$(bindir)' && rm -f" $$files ")"; \
	cd "$(DESTDIR)$(bindir)" && rm -f $$files

clean-binPROGRAMS:
	@list='$(bin_PROGRAMS)'; test -n "$$list" || exit 0; \
	echo " rm -f" $$list; \
	rm -f $$list || exit $$?; \
	test -n "$(EXEEXT)" || exit 0; \
	list=`for p in $$list; do echo "$$p"; done | sed 's/$(EXEEXT)$$//'`; \
	echo " rm -f" $$list; \
	rm -f $$list
atompaw$(EXEEXT): $(atompaw_OBJECTS) $(atompaw_DEPENDENCIES) $(EXTRA_atompaw_DEPENDENCIES) 
	@rm -f atompaw$(EXEEXT)
	$(FCLINK) $(atompaw_OBJECTS) $(atompaw_LDADD) $(LIBS)
graphatom$(EXEEXT): $(graphatom_OBJECTS) $(graphatom_DEPENDENCIES) $(EXTRA_graphatom_DEPENDENCIES) 
	@rm -f graphatom$(EXEEXT)
	$(FCLINK) $(graphatom_OBJECTS) $(graphatom_LDADD) $(LIBS)

mostlyclean-compile:
	-rm -f *.$(OBJEXT)

distclean-compile:
	-rm -f *.tab.c

.F90.o:
	$(PPFCCOMPILE) -c -o $@ $<

.F90.obj:
	$(PPFCCOMPILE) -c -o $@ `$(CYGPATH_W) '$<'`

.F90.lo:
	$(LTPPFCCOMPILE) -c -o $@ $<

.f90.o:
	$(FCCOMPILE) -c -o $@ $<

.f90.obj:
	$(FCCOMPILE) -c -o $@ `$(CYGPATH_W) '$<'`

.f90.lo:
	$(LTFCCOMPILE) -c -o $@ $<

mostlyclean-libtool:
	-rm -f *.lo

clean-libtool:
	-rm -rf .libs _libs

ID: $(HEADERS) $(SOURCES) $(LISP) $(TAGS_FILES)
	list='$(SOURCES) $(HEADERS) $(LISP) $(TAGS_FILES)'; \
	unique=`for i in $$list; do \
	    if test -f "$$i"; then echo $$i; else echo $(srcdir)/$$i; fi; \
	  done | \
	  $(AWK) '{ files[$$0] = 1; nonempty = 1; } \
	      END { if (nonempty) { for (i in files) print i; }; }'`; \
	mkid -fID $$unique
tags: TAGS

TAGS:  $(HEADERS) $(SOURCES)  $(TAGS_DEPENDENCIES) \
		$(TAGS_FILES) $(LISP)
	set x; \
	here=`pwd`; \
	list='$(SOURCES) $(HEADERS)  $(LISP) $(TAGS_FILES)'; \
	unique=`for i in $$list; do \
	    if test -f "$$i"; then echo $$i; else echo $(srcdir)/$$i; fi; \
	  done | \
	  $(AWK) '{ files[$$0] = 1; nonempty = 1; } \
	      END { if (nonempty) { for (i in files) print i; }; }'`; \
	shift; \
	if test -z "$(ETAGS_ARGS)$$*$$unique"; then :; else \
	  test -n "$$unique" || unique=$$empty_fix; \
	  if test $$# -gt 0; then \
	    $(ETAGS) $(ETAGSFLAGS) $(AM_ETAGSFLAGS) $(ETAGS_ARGS) \
	      "$$@" $$unique; \
	  else \
	    $(ETAGS) $(ETAGSFLAGS) $(AM_ETAGSFLAGS) $(ETAGS_ARGS) \
	      $$unique; \
	  fi; \
	fi
ctags: CTAGS
CTAGS:  $(HEADERS) $(SOURCES)  $(TAGS_DEPENDENCIES) \
		$(TAGS_FILES) $(LISP)
	list='$(SOURCES) $(HEADERS)  $(LISP) $(TAGS_FILES)'; \
	unique=`for i in $$list; do \
	    if test -f "$$i"; then echo $$i; else echo $(srcdir)/$$i; fi; \
	  done | \
	  $(AWK) '{ files[$$0] = 1; nonempty = 1; } \
	      END { if (nonempty) { for (i in files) print i; }; }'`; \
	test -z "$(CTAGS_ARGS)$$unique" \
	  || $(CTAGS) $(CTAGSFLAGS) $(AM_CTAGSFLAGS) $(CTAGS_ARGS) \
	     $$unique

GTAGS:
	here=`$(am__cd) $(top_builddir) && pwd` \
	  && $(am__cd) $(top_srcdir) \
	  && gtags -i $(GTAGS_ARGS) "$$here"

cscopelist:  $(HEADERS) $(SOURCES) $(LISP)
	list='$(SOURCES) $(HEADERS) $(LISP)'; \
	case "$(srcdir)" in \
	  [\\/]* | ?:[\\/]*) sdir="$(srcdir)" ;; \
	  *) sdir=$(subdir)/$(srcdir) ;; \
	esac; \
	for i in $$list; do \
	  if test -f "$$i"; then \
	    echo "$(subdir)/$$i"; \
	  else \
	    echo "$$sdir/$$i"; \
	  fi; \
	done >> $(top_builddir)/cscope.files

distclean-tags:
	-rm -f TAGS ID GTAGS GRTAGS GSYMS GPATH tags

distdir: $(DISTFILES)
	@srcdirstrip=`echo "$(srcdir)" | sed 's/[].[^$$\\*]/\\\\&/g'`; \
	topsrcdirstrip=`echo "$(top_srcdir)" | sed 's/[].[^$$\\*]/\\\\&/g'`; \
	list='$(DISTFILES)'; \
	  dist_files=`for file in $$list; do echo $$file; done | \
	  sed -e "s|^$$srcdirstrip/||;t" \
	      -e "s|^$$topsrcdirstrip/|$(top_builddir)/|;t"`; \
	case $$dist_files in \
	  */*) $(MKDIR_P) `echo "$$dist_files" | \
			   sed '/\//!d;s|^|$(distdir)/|;s,/[^/]*$$,,' | \
			   sort -u` ;; \
	esac; \
	for file in $$dist_files; do \
	  if test -f $$file || test -d $$file; then d=.; else d=$(srcdir); fi; \
	  if test -d $$d/$$file; then \
	    dir=`echo "/$$file" | sed -e 's,/[^/]*$$,,'`; \
	    if test -d "$(distdir)/$$file"; then \
	      find "$(distdir)/$$file" -type d ! -perm -700 -exec chmod u+rwx {} \;; \
	    fi; \
	    if test -d $(srcdir)/$$file && test $$d != $(srcdir); then \
	      cp -fpR $(srcdir)/$$file "$(distdir)$$dir" || exit 1; \
	      find "$(distdir)/$$file" -type d ! -perm -700 -exec chmod u+rwx {} \;; \
	    fi; \
	    cp -fpR $$d/$$file "$(distdir)$$dir" || exit 1; \
	  else \
	    test -f "$(distdir)/$$file" \
	    || cp -p $$d/$$file "$(distdir)/$$file" \
	    || exit 1; \
	  fi; \
	done
check-am: all-am
check: check-am
all-am: Makefile $(LTLIBRARIES) $(PROGRAMS)
install-binPROGRAMS: install-libLTLIBRARIES

installdirs:
	for dir in "$(DESTDIR)$(libdir)" "$(DESTDIR)$(bindir)"; do \
	  test -z "$$dir" || $(MKDIR_P) "$$dir"; \
	done
install: install-am
install-exec: install-exec-am
install-data: install-data-am
uninstall: uninstall-am

install-am: all-am
	@$(MAKE) $(AM_MAKEFLAGS) install-exec-am install-data-am

installcheck: installcheck-am
install-strip:
	if test -z '$(STRIP)'; then \
	  $(MAKE) $(AM_MAKEFLAGS) INSTALL_PROGRAM="$(INSTALL_STRIP_PROGRAM)" \
	    install_sh_PROGRAM="$(INSTALL_STRIP_PROGRAM)" INSTALL_STRIP_FLAG=-s \
	      install; \
	else \
	  $(MAKE) $(AM_MAKEFLAGS) INSTALL_PROGRAM="$(INSTALL_STRIP_PROGRAM)" \
	    install_sh_PROGRAM="$(INSTALL_STRIP_PROGRAM)" INSTALL_STRIP_FLAG=-s \
	    "INSTALL_PROGRAM_ENV=STRIPPROG='$(STRIP)'" install; \
	fi
mostlyclean-generic:

clean-generic:

distclean-generic:
	-test -z "$(CONFIG_CLEAN_FILES)" || rm -f $(CONFIG_CLEAN_FILES)
	-test . = "$(srcdir)" || test -z "$(CONFIG_CLEAN_VPATH_FILES)" || rm -f $(CONFIG_CLEAN_VPATH_FILES)

maintainer-clean-generic:
	@echo "This command is intended for maintainers to use"
	@echo "it deletes files that may require special tools to rebuild."
clean: clean-am

clean-am: clean-binPROGRAMS clean-generic clean-libLTLIBRARIES \
	clean-libtool clean-local mostlyclean-am

distclean: distclean-am
	-rm -f Makefile
distclean-am: clean-am distclean-compile distclean-generic \
	distclean-tags

dvi: dvi-am

dvi-am:

html: html-am

html-am:

info: info-am

info-am:

install-data-am:

install-dvi: install-dvi-am

install-dvi-am:

install-exec-am: install-binPROGRAMS install-libLTLIBRARIES

install-html: install-html-am

install-html-am:

install-info: install-info-am

install-info-am:

install-man:

install-pdf: install-pdf-am

install-pdf-am:

install-ps: install-ps-am

install-ps-am:

installcheck-am:

maintainer-clean: maintainer-clean-am
	-rm -f Makefile
maintainer-clean-am: distclean-am maintainer-clean-generic

mostlyclean: mostlyclean-am

mostlyclean-am: mostlyclean-compile mostlyclean-generic \
	mostlyclean-libtool

pdf: pdf-am

pdf-am:

ps: ps-am

ps-am:

uninstall-am: uninstall-binPROGRAMS uninstall-libLTLIBRARIES

.MAKE: install-am install-strip

.PHONY: CTAGS GTAGS all all-am check check-am clean clean-binPROGRAMS \
	clean-generic clean-libLTLIBRARIES clean-libtool clean-local \
	cscopelist ctags distclean distclean-compile distclean-generic \
	distclean-libtool distclean-tags distdir dvi dvi-am html \
	html-am info info-am install install-am install-binPROGRAMS \
	install-data install-data-am install-dvi install-dvi-am \
	install-exec install-exec-am install-html install-html-am \
	install-info install-info-am install-libLTLIBRARIES \
	install-man install-pdf install-pdf-am install-ps \
	install-ps-am install-strip installcheck installcheck-am \
	installdirs maintainer-clean maintainer-clean-generic \
	mostlyclean mostlyclean-compile mostlyclean-generic \
	mostlyclean-libtool pdf pdf-am ps ps-am tags uninstall \
	uninstall-am uninstall-binPROGRAMS uninstall-libLTLIBRARIES


                    # ------------------------------------ #

# Local cleaning
clean-local:
	rm -f *.mod *.MOD *.obj

# Explicit dependencies

abinitinterface.$(LTOBJEXT): atomdata.$(LTOBJEXT) globalmath.$(LTOBJEXT) \
	gridmod.$(LTOBJEXT) excor.$(LTOBJEXT) pseudo.$(LTOBJEXT) \
	interpolation_mod.$(LTOBJEXT) pkginfo.$(LTOBJEXT) libxc_mod.$(LTOBJEXT)

aeatom.$(LTOBJEXT): atomdata.$(LTOBJEXT) excor.$(LTOBJEXT) exx_mod.$(LTOBJEXT) \
	general_mod.$(LTOBJEXT) globalmath.$(LTOBJEXT) gridmod.$(LTOBJEXT) \
	hf_mod.$(LTOBJEXT) ldagga_mod.$(LTOBJEXT)

anderson_driver.$(LTOBJEXT): globalmath.$(LTOBJEXT)

atomdata.$(LTOBJEXT):

atompaw_lib.$(LTOBJEXT):

atompaw_prog.$(LTOBJEXT): globalmath.$(LTOBJEXT) aeatom.$(LTOBJEXT) \
	atomdata.$(LTOBJEXT) atompaw_report.$(LTOBJEXT)  pseudo.$(LTOBJEXT) \
        ini.$(LTOBJEXT) optreal.(LTOBJEXT) radial.$(LTOBJEXT) setlocalpp.$(LTOBJEXT) \
	pseudo_vasp.$(LTOBJEXT) cl_shift.$(LTOBJEXT) rhfatm.$(LTOBJEXT) vasp_pseudo.$(LTOBJEXT) \
	abinitinterface.$(LTOBJEXT) pwscfinterface.$(LTOBJEXT) xmlinterface.$(LTOBJEXT) \
	pkginfo.$(LTOBJEXT) libxc_mod.$(LTOBJEXT)

atompaw_report.$(LTOBJEXT): atomdata.$(LTOBJEXT) fock.$(LTOBJEXT) \
	globalmath.$(LTOBJEXT) gridmod.$(LTOBJEXT) pseudo.$(LTOBJEXT) \
	pseudodata.$(LTOBJEXT) libxc_mod.$(LTOBJEXT) #vasp_pseudo.$(LTOBJEXT)

blockdavidson_mod.$(LTOBJEXT): globalmath.$(LTOBJEXT) search_sort.$(LTOBJEXT)

excor.$(LTOBJEXT): atomdata.$(LTOBJEXT) globalmath.$(LTOBJEXT) gridmod.$(LTOBJEXT) \
	libxc_mod.$(LTOBJEXT)

exx_mod.$(LTOBJEXT): anderson_driver.$(LTOBJEXT) atomdata.$(LTOBJEXT) \
	exxdata.$(LTOBJEXT) fock.$(LTOBJEXT) general_mod.$(LTOBJEXT) \
	globalmath.$(LTOBJEXT) gridmod.$(LTOBJEXT) hf_mod.$(LTOBJEXT) report_mod.$(LTOBJEXT)

exx_pseudo.$(LTOBJEXT): anderson_driver.$(LTOBJEXT) atomdata.$(LTOBJEXT) \
	exx_mod.$(LTOBJEXT) fock.$(LTOBJEXT) globalmath.$(LTOBJEXT) \
	gridmod.$(LTOBJEXT) pseudodata.$(LTOBJEXT) pseudo_sub.$(LTOBJEXT)

exxdata.$(LTOBJEXT):

fock.$(LTOBJEXT): atomdata.$(LTOBJEXT) gridmod.$(LTOBJEXT) globalmath.$(LTOBJEXT)

general_mod.$(LTOBJEXT): atomdata.$(LTOBJEXT) numerov_mod.$(LTOBJEXT) \
	globalmath.$(LTOBJEXT) gridmod.$(LTOBJEXT) radialsr.$(LTOBJEXT)

globalmath.$(LTOBJEXT):

graphatom.$(LTOBJEXT): globalmath.$(LTOBJEXT) aeatom.$(LTOBJEXT) \
	atomdata.$(LTOBJEXT) graphatom_report.$(LTOBJEXT)

graphatom_report.$(LTOBJEXT): globalmath.$(LTOBJEXT) atomdata.$(LTOBJEXT) \
	gridmod.$(LTOBJEXT) report_mod.$(LTOBJEXT)

gridmod.$(LTOBJEXT): globalmath.$(LTOBJEXT)

hf_mod.$(LTOBJEXT): anderson_driver.$(LTOBJEXT) atomdata.$(LTOBJEXT) \
	fock.$(LTOBJEXT) general_mod.$(LTOBJEXT) globalmath.$(LTOBJEXT) \
	gridmod.$(LTOBJEXT) report_mod.$(LTOBJEXT)

hf_pseudo.$(LTOBJEXT): atomdata.$(LTOBJEXT) hf_mod.$(LTOBJEXT) fock.$(LTOBJEXT) \
	globalmath.$(LTOBJEXT) gridmod.$(LTOBJEXT) paw_sub.$(LTOBJEXT) \
	pseudodata.$(LTOBJEXT) pseudo_sub.$(LTOBJEXT)

interpolation_mod.$(LTOBJEXT):

ldagga_mod.$(LTOBJEXT): atomdata.$(LTOBJEXT) anderson_driver.$(LTOBJEXT) \
	excor.$(LTOBJEXT) general_mod.$(LTOBJEXT) globalmath.$(LTOBJEXT) \
	gridmod.$(LTOBJEXT) report_mod.$(LTOBJEXT)

libxc_mod.$(LTOBJEXT): globalmath.$(LTOBJEXT)

numerov_mod.$(LTOBJEXT): gridmod.$(LTOBJEXT) blockdavidson_mod.$(LTOBJEXT)

paw_sub.$(LTOBJEXT): globalmath.$(LTOBJEXT) gridmod.$(LTOBJEXT)

pkginfo.$(LTOBJEXT):

pseudo.$(LTOBJEXT): globalmath.$(LTOBJEXT) atomdata.$(LTOBJEXT) aeatom.$(LTOBJEXT) \
	excor.$(LTOBJEXT) exx_pseudo.$(LTOBJEXT) hf_pseudo.$(LTOBJEXT) \
	numerov_mod.$(LTOBJEXT) paw_sub.$(LTOBJEXT) pseudodata.$(LTOBJEXT) \
	pseudo_sub.$(LTOBJEXT) radialsr.$(LTOBJEXT)

pseudo_sub.$(LTOBJEXT): globalmath.$(LTOBJEXT) gridmod.$(LTOBJEXT) \
	pseudodata.$(LTOBJEXT) excor.$(LTOBJEXT) fock.$(LTOBJEXT)

pseudodata.$(LTOBJEXT): gridmod.$(LTOBJEXT) atomdata.$(LTOBJEXT)

pwscfinterface.$(LTOBJEXT): globalmath.$(LTOBJEXT) atomdata.$(LTOBJEXT) \
	excor.$(LTOBJEXT) gridmod.$(LTOBJEXT) interpolation_mod.$(LTOBJEXT) \
	pseudo.$(LTOBJEXT)

radialsr.$(LTOBJEXT): globalmath.$(LTOBJEXT) gridmod.$(LTOBJEXT) atomdata.$(LTOBJEXT)

report_mod.$(LTOBJEXT): atomdata.$(LTOBJEXT) gridmod.$(LTOBJEXT) excor.$(LTOBJEXT)

search_sort.$(LTOBJEXT):

xmlinterface.$(LTOBJEXT): atomdata.$(LTOBJEXT) globalmath.$(LTOBJEXT) \
	gridmod.$(LTOBJEXT) excor.$(LTOBJEXT) pseudo.$(LTOBJEXT) pkginfo.$(LTOBJEXT)

ini.$(LTOBJEXT):

asa.$(LTOBJEXT): ini.$(LTOBJEXT)

optreal.$(LTOBJEXT): ini.$(LTOBJEXT)

xclib.$(LTOBJEXT): ini.$(LTOBJEXT) 

xclib_grad.$(LTOBJEXT): ini.$(LTOBJEXT) xclib.$(LTOBJEXT)

radial.$(LTOBJEXT): ini.$(LTOBJEXT) xclib.$(LTOBJEXT)

setlocalpp.$(LTOBJEXT): ini.$(LTOBJEXT) radial.$(LTOBJEXT)

pseudo_vasp.$(LTOBJEXT): ini.$(LTOBJEXT) radial.$(LTOBJEXT) setlocalpp.$(LTOBJEXT) \
	optreal.$(LTOBJEXT)

cl_shift.$(LTOBJEXT): pseudo_vasp.$(LTOBJEXT) ini.$(LTOBJEXT) radial.$(LTOBJEXT)

rhfatm.$(LTOBJEXT): pseudo_vasp.$(LTOBJEXT) ini.$(LTOBJEXT) radial.$(LTOBJEXT)  \
	cl_shift.$(LTOBJEXT) xclib.$(LTOBJEXT) 

vasp_pseudo.$(LTOBJEXT): globalmath.$(LTOBJEXT) atomdata.$(LTOBJEXT) aeatom.$(LTOBJEXT) \
	excor.$(LTOBJEXT) exx_pseudo.$(LTOBJEXT) hf_pseudo.$(LTOBJEXT) \
	numerov_mod.$(LTOBJEXT) paw_sub.$(LTOBJEXT) pseudodata.$(LTOBJEXT) \
	pseudo_sub.$(LTOBJEXT) radialsr.$(LTOBJEXT) pseudo.$(LTOBJEXT) \
	ini.$(LTOBJEXT) asa.$(LTOBJEXT) xclib.$(LTOBJEXT) radial.$(LTOBJEXT)  \
	setlocalpp.$(LTOBJEXT) pseudo_vasp.$(LTOBJEXT) \
	atompaw_report.$(LTOBJEXT)


# Tell versions [3.59,3.63) of GNU make to not export all variables.
# Otherwise a system limit (for SysV at least) may be exceeded.
.NOEXPORT:
