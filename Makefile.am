ACLOCAL_AMFLAGS = -I m4

nodist_bin_SCRIPTS = \
  pigx_bs

dist_pkglibexec_SCRIPTS =			\
  BSseq_pipeline.py

pkglibexec_scriptsdir = $(pkglibexecdir)/scripts
dist_pkglibexec_scripts_SCRIPTS =	\
  scripts/fetch_refGene.R			\
  scripts/func_defs.py				\
  scripts/ideoDMC.R					\
  scripts/integrate2finalreport.R	\
  scripts/multireport_functions.R	\
  scripts/report_functions.R

pkglibexec_rulesdir = $(pkglibexecdir)/rules
dist_pkglibexec_rules_SCRIPTS =		\
  rules/post_mapping.rules			\
  rules/samtools.rules				\
  rules/SRA2fastq

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
  etc/programs.json

dist_doc_DATA = 					\
  README.md

CLEANFILES = $(nodist_bin_SCRIPTS)

# This is a hack to embed the full path to the logo after installation.
report_templates/pigx_bsseq_logo.html: report_templates/pigx_bsseq_logo.html.in
	sed -e "s|#:MAKE|$(pkgdatadir)|" $< > $@
