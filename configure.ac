dnl -*- Autoconf -*-

AC_INIT([PiGx_rnaseq], [m4_translit(m4_esyscmd([cat VERSION]),m4_newline)])
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
find_or_override_prog([TRIMGALORE],      [trim_galore])
find_or_override_prog([STAR],            [STAR])
find_or_override_prog([MULTIQC],         [multiqc])
find_or_override_prog([HTSEQ_COUNT],     [htseq-count])
find_or_override_prog([BAMCOVERAGE],     [bamCoverage]) dnl deeptools

find_or_override_prog([SALMON],     [salmon])

find_or_override_prog([R],               [R])
find_or_override_prog([RSCRIPT],         [Rscript])


AC_ARG_ENABLE([r-packages-check],
  AS_HELP_STRING([--disable-r-packages-check], [Do not check any R packages.]),
  [AC_MSG_NOTICE([Skipping R packages check.  Be careful!])],
  [dnl
AX_R_PACKAGE([rmarkdown])
AX_R_PACKAGE([knitr])
AX_R_PACKAGE([ggplot2])
AX_R_PACKAGE([ggrepel])
AX_R_PACKAGE([DESeq2])
AX_R_PACKAGE([DT])
AX_R_PACKAGE([pheatmap])
AX_R_PACKAGE([corrplot])
AX_R_PACKAGE([reshape2])
AX_R_PACKAGE([plotly])
AX_R_PACKAGE([scales])
AX_R_PACKAGE([crosstalk])
AX_R_PACKAGE([gProfileR])
AX_R_PACKAGE([rtracklayer])
AX_R_PACKAGE([SummarizedExperiment])
AX_R_PACKAGE([tximport])
AX_R_PACKAGE([rjson])
  ])


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

AC_CONFIG_FILES([etc/settings.yaml
                 Makefile])
AC_CONFIG_FILES([test.sh], [chmod +x test.sh])
AC_CONFIG_FILES([tests/test_trim_galore/test.sh], [chmod +x tests/test_trim_galore/test.sh])
AC_CONFIG_FILES([tests/test_star_map/test.sh], [chmod +x tests/test_star_map/test.sh])
AC_CONFIG_FILES([pigx-rnaseq], [chmod +x pigx-rnaseq])

AC_OUTPUT
