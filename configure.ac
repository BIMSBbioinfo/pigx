dnl -*- Autoconf -*-

AC_INIT([PiGx_scrnaseq], [m4_translit(m4_esyscmd([cat VERSION]),m4_newline)])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([gnu color-tests tar-pax -Wall -Wno-portability foreign])
AM_SILENT_RULES([yes])

dnl We need sed for version checks
AC_PROG_SED

dnl TODO: this is not pretty.  We shouldn't need to do this as long as
dnl these jars are on the CLASSPATH.

dnl Load Python support for Automake and ensure we have a Python
dnl interpreter for at least version 3.6.  This sets the "PYTHON"
dnl variable.
dnl https://www.gnu.org/software/automake/manual/html_node/Python.html
AM_PATH_PYTHON([3.6])

AX_PYTHON_MODULE([yaml],   "required")
AX_PYTHON_MODULE([pandas], "required")
AX_PYTHON_MODULE([numpy],  "required")
AX_PYTHON_MODULE([loompy], "required")
AX_PYTHON_MODULE([magic], "required")

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

find_or_override_prog([GNUBASH],         [bash])
find_or_override_prog([SNAKEMAKE],       [snakemake])
find_or_override_prog([PANDOC],          [pandoc])
find_or_override_prog([PANDOC_CITEPROC], [pandoc-citeproc])
find_or_override_prog([FASTQC],          [fastqc])
find_or_override_prog([STAR],            [STAR])
find_or_override_prog([SAMTOOLS],        [samtools])
find_or_override_prog([LN],              [ln])
find_or_override_prog([CAT],             [cat])
find_or_override_prog([ZCAT],            [zcat])
find_or_override_prog([RM],              [rm])
find_or_override_prog([TOUCH],           [touch])
find_or_override_prog([PERL],            [perl])
find_or_override_prog([R],               [R])
find_or_override_prog([RSCRIPT],         [Rscript])
find_or_override_prog([JAVA],            [java])
find_or_override_prog([FLEXBAR],         [flexbar])
find_or_override_prog([JELLYFISH],       [jellyfish])

AC_ARG_ENABLE([environment-capture],
  AS_HELP_STRING([--disable-environment-capture], [Do not capture environment variables.]),
  [AC_MSG_NOTICE([Environment variables will not be captured.])],
  [dnl
AC_MSG_NOTICE([Environment variables will be captured.])
capture_environment=yes
AC_ARG_VAR(R_LIBS_SITE, The search path for all R packages)
AC_ARG_VAR(PYTHONPATH, The search path for all Python packages)
  ])

AC_SUBST([capture_environment])

AC_ARG_ENABLE([r-packages-check],
  AS_HELP_STRING([--disable-r-packages-check], [Do not check any R packages.]),
  [AC_MSG_NOTICE([Skipping R packages check.  Be careful!])],
  [dnl
AX_R_PACKAGE([argparser])
AX_R_PACKAGE([dplyr])
AX_R_PACKAGE([dropbead])
AX_R_PACKAGE([GenomicAlignments])
AX_R_PACKAGE([GenomicFiles])
AX_R_PACKAGE([GenomicRanges])
AX_R_PACKAGE([ggplot2])
AX_R_PACKAGE([pheatmap])
AX_R_PACKAGE([Rsamtools])
AX_R_PACKAGE([rtracklayer])
AX_R_PACKAGE([stringr])
AX_R_PACKAGE([yaml])
AX_R_PACKAGE([cowplot])
AX_R_PACKAGE([data.table])
AX_R_PACKAGE([DelayedArray])
AX_R_PACKAGE([DelayedMatrixStats])
AX_R_PACKAGE([DT])
AX_R_PACKAGE([HDF5Array])
AX_R_PACKAGE([rmarkdown])
AX_R_PACKAGE([Rtsne])
AX_R_PACKAGE([scater])
AX_R_PACKAGE([scran])
AX_R_PACKAGE([Seurat])
AX_R_PACKAGE([SingleCellExperiment])
  ])

AC_CONFIG_FILES([etc/settings.yaml
                 scripts/scrnaReport.Rmd
                 Makefile
                 qsub-template.sh])
AC_CONFIG_FILES([tests/test.sh], [chmod +x tests/test.sh])
AC_CONFIG_FILES([pigx-scrnaseq], [chmod +x pigx-scrnaseq])
AC_OUTPUT
