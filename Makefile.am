ACLOCAL_AMFLAGS = -I m4

nodist_bin_SCRIPTS = \
  pigx-rnaseq

dist_pkglibexec_SCRIPTS = \
  pigx_rnaseq.py

pkglibexec_scriptsdir = $(pkglibexecdir)/scripts
dist_pkglibexec_scripts_SCRIPTS =					\
  scripts/counts_matrix_from_STAR.R					\
  scripts/runDeseqReport.R							\
  scripts/translate_sample_sheet_for_report.R		\
  scripts/deseqReport.Rmd		\
  scripts/counts_matrix_from_SALMON.R

dist_pkgdata_DATA =									\
  etc/settings.yaml									\
  etc/pretty.txt									\
  images/Logo_PiGx.png

dist_doc_DATA = 									\
  README.md

EXTRA_DIST =										\
  tests/sample_sheet.csv							\
  tests/sample_data/sample.gtf						\
  tests/sample_data/reads/HBR_Rep1.read2.fastq.gz	\
  tests/sample_data/reads/HBR_Rep1.read1.fastq.gz	\
  tests/sample_data/reads/UHR_Rep1.read1.fastq.gz	\
  tests/sample_data/reads/UHR_Rep3.read1.fastq.gz	\
  tests/sample_data/reads/HBR_Rep2.read1.fastq.gz	\
  tests/sample_data/reads/UHR_Rep2.read1.fastq.gz	\
  tests/sample_data/reads/UHR_Rep2.read2.fastq.gz	\
  tests/sample_data/reads/HBR_Rep3.read2.fastq.gz	\
  tests/sample_data/reads/UHR_Rep1.read2.fastq.gz	\
  tests/sample_data/reads/HBR_Rep3.read1.fastq.gz	\
  tests/sample_data/reads/UHR_Rep3.read2.fastq.gz	\
  tests/sample_data/reads/HBR_Rep2.read2.fastq.gz	\
  tests/sample_data/sample.fasta					\
  tests/settings.yaml

AM_TESTS_ENVIRONMENT = srcdir="$(abs_top_srcdir)"

TESTS = \
  tests/test_salmon/test_salmon_index.sh				\
  tests/test_salmon/test_salmon_quant.sh				\
  tests/test_star/test.sh						\
  tests/test_trim_galore/test.sh					\
  test.sh 

clean-local:
	-rm -rf .snakemake

CLEANFILES = $(nodist_bin_SCRIPTS) config.json

# Build docker image with Guix
docker: dist
	guix pack -e '(load "guix.scm")' --with-source=pigx_rnaseq-$(VERSION).tar.gz -f docker -S /bin=bin -S /lib=lib -S /share=share glibc-utf8-locales tzdata coreutils
