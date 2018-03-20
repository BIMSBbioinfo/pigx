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
  scripts/ChIPQC.R						\
  scripts/ChIPQC_Functions.R			\
  scripts/Check_Config.py				\
  scripts/ChipQC.rmd					\
  scripts/ConstructGenomicWindows.R		\
  scripts/Extract_Nucleotide_Frequency.R\
  scripts/Extract_Signal_Annotation.R	\
  scripts/Extract_Signal_Peaks.R		\
  scripts/Feature_Combination.R			\
  scripts/Functions_Helper.R			\
  scripts/Functions_Sample_Report.R		\
  scripts/Knit_Report.R					\
  scripts/Make_UCSC_HUB.R				\
  scripts/Parse_Bowtie2log.R			\
  scripts/Peak_Statistics.R				\
  scripts/Prepare_Annotation.R			\
  scripts/Sample_Report.rmd				\
  scripts/SnakeFunctions.py				\
  scripts/Summarize_Data_For_Report.R

pkglibexec_rulesdir = $(pkglibexecdir)/Rules

dist_pkglibexec_rules_SCRIPTS =   		\
  Rules/Annotate_Peaks.py       		\
  Rules/BamToBigWig.py          		\
  Rules/ChIPQC.py           			\
  Rules/Extract_Signal_Annotation.py  	\
  Rules/FastQC.py           			\
  Rules/Feature_Combination.py      	\
  Rules/IDR.py              			\
  Rules/Knit_Report.py           		\
  Rules/Mapping.py            			\
  Rules/MultiQC.py            			\
  Rules/Parse_Bowtie2log.py            	\
  Rules/Peak_Calling.py         		\
  Rules/Peak_Statistics.py         		\
  Rules/Prepare_Annotation.py     		\
  Rules/Summarize_Data_For_Report.py	\
  Rules/Trimming.py                		\
  Rules/UCSC_Hub.py

dist_pkgdata_DATA =						\
  etc/settings.yaml						\
  etc/pretty.txt						\
  etc/sample_sheet.csv.example  

AM_TESTS_ENVIRONMENT = srcdir="$(abs_top_srcdir)" builddir="$(abs_top_builddir)" PIGX_UNINSTALLED=1 PIGX_UGLY=1 SOURCE_DATE_EPOCH=1

TESTS = \
  Tests/test.sh

EXTRA_DIST =												\
  Tests/sample_sheet.csv									\
  Tests/settings.yaml										\
  Tests/genome/genome_subset_mm9_chr1_13590945_17791527.fa.gz	\
  Tests/genome/mm_chr19.gtf									\
  Tests/in/ChIP.fq.gz										\
  Tests/in/ChIP_1.fq.gz										\
  Tests/in/ChIP_2.fq.gz										\
  Tests/in/Cont.fq.gz										\
  Tests/in/ChIPpe_R1.fq.gz									\
  Tests/in/ChIPpe_R2.fq.gz

clean-local:
	-rm -rf .snakemake
	-rm -rf Tests/out/
	-rm -rf Tests/out/.snakemake

CLEANFILES = $(nodist_bin_SCRIPTS) config.json

# Build docker image with Guix
docker: dist
	guix pack -e '(load "guix.scm")' --with-source=pigx_chipseq-$(VERSION).tar.gz -f docker -S /bin=bin -S /lib=lib -S /share=share glibc-utf8-locales tzdata coreutils bash
