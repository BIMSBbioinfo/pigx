ACLOCAL_AMFLAGS = -I m4

nodist_bin_SCRIPTS = \
  pigx_bs

nodist_libexec_SCRIPTS =			\
  scripts/fetch_refGene.R			\
  scripts/func_defs.py				\
  scripts/ideoDMC.R					\
  scripts/integrate2finalreport.R	\
  scripts/multireport_functions.R	\
  scripts/report_functions.R		\
  rules/post_mapping.rules			\
  rules/samtools.rules				\
  rules/SRA2fastq

dist_doc_DATA = 					\
  README.md

CLEANFILES = $(nodist_bin_SCRIPTS)

# This is a hack to embed the full path to the logo after installation.
report_templates/pigx_bsseq_logo.html: report_templates/pigx_bsseq_logo.html.in
	sed -e "s|#:MAKE|$(pkgdatadir)|" $< > $@
