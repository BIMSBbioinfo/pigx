ACLOCAL_AMFLAGS = -I m4

nodist_bin_SCRIPTS = \
  pigx-chipseq

dist_pkglibexec_SCRIPTS = \
  Snake_ChIPseq.py

pkglibexec_scriptsdir = $(pkglibexecdir)/scripts
dist_pkglibexec_scripts_SCRIPTS =		\
  scripts/Annotate_Peaks.R				\
  scripts/Argument_Parser.R				\
  scripts/BigWigExtend.R				\
  scripts/Check_Config.py				\
  scripts/ChIPQC.R						\
  scripts/ChipQC.rmd					\
  scripts/Extract_Signal_Annotation.R	\
  scripts/Extract_Signal_Peaks.R		\
  scripts/Feature_Combination.R			\
  scripts/Functions_Helper.R			\
  scripts/Make_UCSC_HUB.R				\
  scripts/Prepare_Annotation.R			\
  scripts/Sample_Report.rmd				\
  scripts/SnakeFunctions.py

dist_pkgdata_DATA =						\
  etc/settings.yaml						\
  etc/pretty.txt						\
  Rules/Annotate_Peaks.py				\
  Rules/BamToBigWig.py					\
  Rules/ChIPQC.py						\
  Rules/Extract_Signal_Annotation.py	\
  Rules/FastQC.py						\
  Rules/Feature_Combination.py			\
  Rules/IDR.py							\
  Rules/Mapping.py						\
  Rules/Peak_Calling.py					\
  Rules/Prepare_Annotation.py			\
  Rules/UCSC_Hub.py

# Build docker image with Guix
docker: dist
	guix pack -e '(load "guix.scm")' --with-source=pigx_chipseq-$(VERSION).tar.gz -f docker -S /bin=bin -S /lib=lib -S /share=share glibc-utf8-locales tzdata coreutils
