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
  tests/test_deseq_reports/test.sh					\
  tests/test_multiqc/test.sh						\
  tests/test_star_counts/test.sh					\
  tests/test_fastqc/test.sh						\
  tests/test_salmon_counts/test.sh					\
  tests/test_salmon/test_salmon_index.sh				\
  tests/test_salmon/test_salmon_quant.sh				\
  tests/test_star/test.sh						\
  tests/test_trim_galore/test.sh					\
  test.sh 

clean-local:
	-rm -rf .snakemake
	-rm -rf tests/output
	-rm -rf tests/output/.snakemake
	-rm -rf tests/test_deseq_reports/config.json
	-rm -rf tests/test_deseq_reports/logs
	-rm -rf tests/test_deseq_reports/pigx_work
	-rm -rf tests/test_deseq_reports/report
	-rm -rf tests/test_deseq_reports/.snakemake
	-rm -rf tests/test_deseq_reports/test.sh
	-rm -rf tests/test_fastqc/fastqc/
	-rm -rf tests/test_fastqc/logs/
	-rm -rf tests/test_fastqc/pigx_work/
	-rm -rf tests/test_fastqc/.snakemake
	-rm -rf tests/test_fastqc/test.sh
	-rm -rf tests/test_fastqc/test.sh.log
	-rm -rf tests/test_fastqc/test.sh.trs
	-rm -rf tests/test_multiqc/config.json
	-rm -rf tests/test_multiqc/logs
	-rm -rf tests/test_multiqc/multiqc
	-rm -rf tests/test_multiqc/pigx_work
	-rm -rf tests/test_multiqc/.snakemake
	-rm -rf tests/test_multiqc/test.sh
	-rm -rf tests/test_salmon_counts/colData.tsv
	-rm -rf tests/test_salmon_counts/logs/
	-rm -rf tests/test_salmon_counts/pigx_work/
	-rm -rf tests/test_salmon_counts/salmon_output/
	-rm -rf tests/test_salmon_counts/.snakemake
	-rm -rf tests/test_salmon_counts/test.sh
	-rm -rf tests/test_salmon_counts/test.sh.log
	-rm -rf tests/test_salmon_counts/test.sh.trs
	-rm -rf tests/test_salmon/logs/
	-rm -rf tests/test_salmon/pigx_work/
	-rm -rf tests/test_salmon/salmon_index/
	-rm -rf tests/test_salmon/salmon_output/
	-rm -rf tests/test_salmon/.snakemake
	-rm -rf tests/test_salmon/test_salmon_index.sh
	-rm -rf tests/test_salmon/test_salmon_index.sh
	-rm -rf tests/test_salmon/test_salmon_index.sh.log
	-rm -rf tests/test_salmon/test_salmon_index.sh.trs
	-rm -rf tests/test_salmon/test_salmon_quant.sh
	-rm -rf tests/test_salmon/test_salmon_quant.sh
	-rm -rf tests/test_salmon/test_salmon_quant.sh.log
	-rm -rf tests/test_salmon/test_salmon_quant.sh.trs
	-rm -rf tests/test_star_counts/config.json
	-rm -rf tests/test_star_counts/pigx_work
	-rm -rf tests/test_star_counts/preprocessed_data
	-rm -rf tests/test_star_counts/.snakemake
	-rm -rf tests/test_star_counts/test.sh
	-rm -rf tests/test_star/Log.out
	-rm -rf tests/test_star/logs/
	-rm -rf tests/test_star/mapped_reads/
	-rm -rf tests/test_star/pigx_work/
	-rm -rf tests/test_star/.snakemake
	-rm -rf tests/test_star/star_index/
	-rm -rf tests/test_star/test.sh
	-rm -rf tests/test_star/test.sh.log
	-rm -rf tests/test_star/test.sh.trs
	-rm -rf tests/test_trim_galore/output_paired_end/
	-rm -rf tests/test_trim_galore/output_single_end/
	-rm -rf tests/test_trim_galore/test.sh
	-rm -rf tests/test_trim_galore/test.sh.log
	-rm -rf tests/test_trim_galore/test.sh.trs		
CLEANFILES = $(nodist_bin_SCRIPTS) config.json		
# Build docker image with Guix
docker: dist
	guix pack -e '(load "guix.scm")' --with-source=pigx_rnaseq-$(VERSION).tar.gz -f docker -S /bin=bin -S /lib=lib -S /share=share glibc-utf8-locales tzdata coreutils
