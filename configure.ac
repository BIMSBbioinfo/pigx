dnl -*- Autoconf -*-

AC_INIT([PiGx], [m4_translit(m4_esyscmd([cat VERSION]),m4_newline)])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([gnu color-tests tar-pax -Wall -Wno-portability foreign])
AM_SILENT_RULES([yes])

dnl We need sed for version checks
AC_PROG_SED

dnl Load Python support for Automake and ensure we have a Python
dnl interpreter for at least version 3.5.  This sets the "PYTHON"
dnl variable.
dnl https://www.gnu.org/software/automake/manual/html_node/Python.html
AM_PATH_PYTHON([3.5])
AX_PYTHON_MODULE([yaml], "required")

AC_ARG_ENABLE([bsseq],
   AS_HELP_STRING([--disable-bsseq], [Do not build the bisulfite sequencing pipeline.]))
if test "x$enable_bsseq" != "xno"; then
  AC_CONFIG_SUBDIRS([pipelines/bsseq])
fi

AC_ARG_ENABLE([rnaseq],
   AS_HELP_STRING([--disable-rnaseq], [Do not build the RNAseq pipeline.]))
if test "x$enable_rnaseq" != "xno"; then
  AC_CONFIG_SUBDIRS([pipelines/rnaseq])
fi

AC_ARG_ENABLE([chipseq],
   AS_HELP_STRING([--enable-chipseq], [Build the chipseq pipeline.]))
if test "x$enable_chipseq" = "xyes"; then
  AC_CONFIG_SUBDIRS([pipelines/chipseq])
fi

AC_ARG_ENABLE([crispr],
   AS_HELP_STRING([--enable-crispr], [Build the crispr pipeline.]))
if test "x$enable_crispr" = "xyes"; then
  AC_CONFIG_SUBDIRS([pipelines/crispr])
fi

AC_CONFIG_FILES([pigx], [chmod +x pigx])
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
