ACLOCAL_AMFLAGS = -I m4

nodist_bin_SCRIPTS = \
  pigx-sars-cov2-ww

dist_pkglibexec_SCRIPTS = \
  snakefile.py

pkglibexec_scriptsdir = $(pkglibexecdir)/scripts
dist_pkglibexec_scripts_SCRIPTS =

dist_pkgdata_DATA =            \
  etc/pretty.txt

dist_doc_DATA = 									\
  README.md

EXTRA_DIST = META.in

AM_TESTS_ENVIRONMENT = srcdir="$(abs_top_srcdir)" builddir="$(abs_top_builddir)" PIGX_UNINSTALLED=1 PIGX_UGLY=1

TESTS = \
  tests/tests/vep.sh

clean-local:
	-rm -rf .snakemake

CLEANFILES = $(nodist_bin_SCRIPTS) config.json		

# Build docker image with Guix
docker: dist
	guix pack -C none -e '(load "guix.scm")' --with-source=pigx_sars-cov2-ww-$(VERSION).tar.gz -f docker -S /bin=bin -S /lib=lib -S /share=share glibc-utf8-locales tzdata coreutils bash
