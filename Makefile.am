ACLOCAL_AMFLAGS = -I m4

nodist_bin_SCRIPTS = \
  pigx-bsseq

dist_pkglibexec_SCRIPTS =			\
  BSseq_pipeline.py

pkglibexec_scriptsdir = $(pkglibexecdir)/scripts
dist_pkglibexec_scripts_SCRIPTS =	\
  scripts/generate_report.R			\
  scripts/generate_multireport.R	\
  scripts/fetch_procedures.R		\
  scripts/fetch_refGene.R			\
  scripts/func_defs.py				\
  scripts/ideoDMC.R					\
  scripts/integrate2finalreport.R	\
  scripts/report_functions.R		\
  scripts/export_bw.R

pkgdata_report_templatesdir = $(pkgdatadir)/report_templates
dist_pkgdata_report_templates_SCRIPTS = 			\
  report_templates/annotation.report.Rmd			\
  report_templates/annotation.report.diff.meth.Rmd	\
  report_templates/diffmeth.report.Rmd				\
  report_templates/index.Rmd						\
  report_templates/methCall.report.Rmd				\
  report_templates/methseg.report.Rmd				\
  report_templates/pigx_bsseq_logo.html				\
  report_templates/references.Rmd					\
  report_templates/sessioninfo.Rmd					\
  report_templates/reports.bib

dist_pkgdata_DATA =					\
  etc/settings.yaml					\
  etc/pretty.txt

dist_doc_DATA = 					\
  README.md

EXTRA_DIST =										\
  test/genome/myref_hg38_chr1_7291263_7302041.fasta	\
  test/in/PE_1.fq.gz								\
  test/in/PE_2.fq.gz								\
  test/in/SE_techrep1.fq.gz							\
  test/in/SE_techrep2.fq.gz							\
  test/samplesheet.csv								\
  test/settings.yaml

TESTS = test.sh

clean-local:
	-rm -rf .snakemake

CLEANFILES = $(nodist_bin_SCRIPTS) report_templates/pigx_bsseq_logo.html config.json

# This is a hack to embed the logo in the html file.
report_templates/pigx_bsseq_logo.html:
	@sed "s|#:MAKE|$(shell base64 -w0 < images/Logo_PiGx.png)|" report_templates/pigx_bsseq_logo.html.in > $@

# Build docker image with Guix
docker: dist
	guix pack -e '(load "guix.scm")' --with-source=pigx_bsseq-$(VERSION).tar.gz -f docker -S /bin=bin -S /lib=lib -S /share=share glibc-utf8-locales tzdata coreutils

# For developers: fetch and convert bed files for hg19 assembly.  Requires bedops.
fetch-cpg:
	wget -qO- http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/cpgIslandExt.txt.gz \
	   | gunzip -c \
	   | awk 'BEGIN{ OFS="\t"; }{ print $$2, $$3, $$4, $$5$$6, substr($$0, index($$0, $$7)); }' \
	   | sort-bed - | gzip > cpgIslandExt.hg19.bed.gz

fetch-ref:
	wget -qO- http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/refGene.txt.gz \
	   | gunzip -c \
	   | awk 'BEGIN{ OFS="\t"; }{ print $$3, $$5, $$6, $$2, $$12, $$4, $$7, $$8; }' \
	   | sort-bed - | gzip > refGene.hg19.bed.gz
