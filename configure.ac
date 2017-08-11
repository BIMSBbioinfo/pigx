dnl -*- Autoconf -*-

AC_INIT([PIGx-bsseq], [0.0.1])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
dnl AM_INIT_AUTOMAKE([gnu color-tests tar-pax -Wall -Wno-portability foreign])
AM_SILENT_RULES([yes])

dnl We need sed for version checks
AC_PROG_SED

dnl Check for required programmes and store their full path in the
dnl given variables.  The variables are used to substitute
dnl placeholders in the scripts.

dnl The variables only need to be set manually when they cannot be
dnl found via the PATH environment variable.
AC_DEFUN([find_or_override_prog],
[AC_ARG_VAR($1, override location of $2 executable)dnl
AS_IF([test -z "$$1"],dnl
      [AC_PATH_PROG([$1], [$2], no)dnl
      if test x"$$1" == x"no" ; then
        AC_MSG_ERROR([Could not find $2.])
      fi
      ],dnl
      [AC_MSG_NOTICE([Using $$1 as $2 executable.])])])

dnl Ensure that an R package can be loaded, unless the HURRY
dnl environment variable is set.
AC_DEFUN([AX_TEST_R_PACKAGE], [
    pushdef([PKG],$1)
    AC_MSG_CHECKING([for R package "PKG"])
    if test -z "$HURRY"
    then
      TEST=$( R --silent --vanilla -e 'if(!require("PKG")) stop("not found")' 2>/dev/null )
      AS_IF([test $? != 0],
            [AC_MSG_ERROR([R package PKG could not be loaded])],
            [AC_MSG_RESULT([OK])])
    else
      AC_MSG_RESULT([OK])
    fi
    popdef([PKG])
])

find_or_override_prog([GNUBASH],         [bash])
find_or_override_prog([PANDOC],          [pandoc])
find_or_override_prog([PANDOC_CITEPROC], [pandoc-citeproc])
find_or_override_prog([R],               [R])
find_or_override_prog([RSCRIPT],         [Rscript])

AX_TEST_R_PACKAGE([genomation])
AX_TEST_R_PACKAGE([methylKit])
AX_TEST_R_PACKAGE([AnnotationHub])
AX_TEST_R_PACKAGE([DT])
AX_TEST_R_PACKAGE([rtracklayer])
AX_TEST_R_PACKAGE([rmarkdown])

find_or_override_prog([SAMTOOLS],      [samtools])
AC_ARG_ENABLE([samtools-version-check],
  AS_HELP_STRING([--disable-samtools-version-check], [Do not check the Samtools version.]),
  [AC_MSG_NOTICE([Skipping Samtools version check.  Be careful!])],
  [dnl
SAMTOOLS_VERSION_MIN=1.3
SAMTOOLS_VERSION=$([$SAMTOOLS --version | $SED 's/^samtools //;q'])
AX_COMPARE_VERSION([$SAMTOOLS_VERSION], [lt], [$SAMTOOLS_VERSION_MIN],
      [AC_MSG_ERROR([Please install Samtools $SAMTOOLS_VERSION_MIN.])])
  ])

dnl AC_CONFIG_FILES([Makefile
dnl                  fp_from_trained_params.R
dnl                  generate_cut_matrix.sh.in
dnl                  get_cut_matrix.R
dnl                  make_1bp_bam_for_fp.sh])

dnl AC_CONFIG_FILES([test.sh], [chmod +x test.sh])
dnl AC_CONFIG_FILES([find_footprints_from_trained_parameters.sh.in],
dnl                 [chmod +x find_footprints_from_trained_parameters.sh.in])

AC_OUTPUT
