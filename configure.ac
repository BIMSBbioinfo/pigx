dnl -*- Autoconf -*-

AC_INIT([PiGx], [m4_translit(m4_esyscmd([cat VERSION]),m4_newline)])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([gnu color-tests tar-pax -Wall -Wno-portability foreign])
AM_SILENT_RULES([yes])

dnl We need sed for version checks
AC_PROG_SED

dnl Load Python support for Automake and ensure we have a Python
dnl interpreter for at least version 3.6.  This sets the "PYTHON"
dnl variable.
dnl https://www.gnu.org/software/automake/manual/html_node/Python.html
AM_PATH_PYTHON([3.6])

pipelines=""

AC_ARG_ENABLE([bundle],
   AS_HELP_STRING([--enable-bundle], [Build pipelines in sub-directories.]))

AC_ARG_ENABLE([bsseq],
   AS_HELP_STRING([--disable-bsseq], [Do not build the bisulfite sequencing pipeline.]))
if test "x$enable_bsseq" != "xno"; then
  if test "x$enable_bundle" = "xyes"; then
    AC_CONFIG_SUBDIRS([pipelines/bsseq])
  else
    AC_PATH_PROG([PIGX_BSSEQ], [pigx-bsseq])
  fi
  pipelines="$pipelines bsseq"
fi

AC_ARG_ENABLE([rnaseq],
   AS_HELP_STRING([--disable-rnaseq], [Do not build the RNAseq pipeline.]))
if test "x$enable_rnaseq" != "xno"; then
  if test "x$enable_bundle" = "xyes"; then
    AC_CONFIG_SUBDIRS([pipelines/rnaseq])
  else
    AC_PATH_PROG([PIGX_RNASEQ], [pigx-rnaseq])
  fi
  pipelines="$pipelines rnaseq"
fi

AC_ARG_ENABLE([chipseq],
   AS_HELP_STRING([--disable-chipseq], [Do not build the chipseq pipeline.]))
if test "x$enable_chipseq" != "xno"; then
  if test "x$enable_bundle" = "xyes"; then
    AC_CONFIG_SUBDIRS([pipelines/chipseq])
  else
    AC_PATH_PROG([PIGX_CHIPSEQ], [pigx-chipseq])
  fi
  pipelines="$pipelines chipseq"
fi

AC_ARG_ENABLE([crispr],
   AS_HELP_STRING([--enable-crispr], [Build the crispr pipeline.]))
if test "x$enable_crispr" = "xyes"; then
  if test "x$enable_bundle" = "xyes"; then
    AC_CONFIG_SUBDIRS([pipelines/crispr])
  else
    AC_PATH_PROG([PIGX_CRISPR], [pigx-crispr])
  fi
  pipelines="$pipelines crispr"
fi

AC_ARG_ENABLE([scrnaseq],
   AS_HELP_STRING([--disable-scrnaseq], [Do not build the single cell RNAseq pipeline.]))
if test "x$enable_scrnaseq" != "xno"; then
  if test "x$enable_bundle" = "xyes"; then
    AC_CONFIG_SUBDIRS([pipelines/scrnaseq])
  else
    AC_PATH_PROG([PIGX_SCRNASEQ], [pigx-scrnaseq])
  fi
  pipelines="$pipelines scrnaseq"
fi

AC_SUBST([pipelines], [$pipelines])
AC_CONFIG_FILES([pigx], [chmod +x pigx])
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
