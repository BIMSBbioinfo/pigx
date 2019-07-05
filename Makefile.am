ACLOCAL_AMFLAGS = -I m4

nodist_bin_SCRIPTS = \
  pigx-scrnaseq

dist_pkglibexec_SCRIPTS =			\
  Snake_Dropseq.py \
  qsub-template.sh

pkglibexec_scriptsdir = $(pkglibexecdir)/scripts
dist_pkglibexec_scripts_SCRIPTS =					\
  scripts/Accessory_Functions.py					\
  scripts/Argument_Parser.R						\
  scripts/BamToBigWig.R								\
  scripts/Extract_Downstream_Statistics.R				\
  scripts/Extract_Read_Statistics.R					\
  scripts/Find_Absolute_Read_Cutoff.R				\
  scripts/change_gtf_id.R							\
  scripts/combine_UMI_matrices.py					\
  scripts/convert_loom_to_singleCellExperiment.R	\
  scripts/convert_matrix_to_loom.py					\
  scripts/renderReport.R							\
  scripts/Run_Rscript.py						\
  scripts/Sample_Sheet_Class.py						\
  scripts/scrnaReport.Rmd							\
  scripts/validate_input.py \
  scripts/convert_loom_to_Seurat.R \
  scripts/loom_Functions.R

dist_pkgdata_DATA =					\
  etc/sample_sheet.csv.example		\
  etc/settings.yaml					\
  etc/pretty.txt

dist_doc_DATA = 					\
  README.md

EXTRA_DIST =										   \
  tests/sample_sheet.csv                               \
  tests/settings.yaml                                  \
  tests/sample_data/hg19_chr19.10094960.25108012.fa    \
  tests/sample_data/hg19_chr19.10094960.25108012.gtf   \
  tests/sample_data/reads/HEK_4h_br1_R1.fastq.gz	\
  tests/sample_data/reads/HEK_4h_br1_R2.fastq.gz   \
  tests/sample_data/reads/HEK_0h_br1_R1_1.fastq.gz \
  tests/sample_data/reads/HEK_0h_br1_R2_1.fastq.gz \
  tests/sample_data/reads/HEK_0h_br1_R1_2.fastq.gz \
  tests/sample_data/reads/HEK_0h_br1_R2_2.fastq.gz


AM_TESTS_ENVIRONMENT = srcdir="$(abs_top_srcdir)" builddir="$(abs_top_builddir)" PIGX_UNINSTALLED=1 PIGX_UGLY=1

TESTS = \
  tests/test.sh

EXTRA_DIST += $(TESTS)

clean-local:
	-rm -rf .snakemake

CLEANFILES = $(nodist_bin_SCRIPTS) config.json

# Build docker image with Guix
docker: dist
	guix pack -e '(load "guix.scm")' --with-source=pigx_scrnaseq-$(VERSION).tar.gz -f docker -S /bin=bin -S /lib=lib -S /share=share glibc-utf8-locales tzdata coreutils
