ACLOCAL_AMFLAGS = -I m4

nodist_bin_SCRIPTS = \
  pigx-sars-cov2-ww

dist_pkglibexec_SCRIPTS = \
  snakefile.py

pkglibexec_scriptsdir = $(pkglibexecdir)/scripts
dist_pkglibexec_scripts_SCRIPTS = \
  scripts/download_databases.sh \
  scripts/report_scripts/qc_report_per_sample.Rmd \
  scripts/report_scripts/taxonomic_classification.Rmd \
  scripts/report_scripts/index.Rmd \
  scripts/report_scripts/variantreport_p_sample.Rmd \
  scripts/generateNavigation.R \
  scripts/renderReport.R \
  scripts/get_qc_table.py \
  scripts/parse_vep.py \
  scripts/vcfTocsv.py \
  scripts/vep_rule.py

dist_pkgdata_DATA =            \
  etc/pretty.txt \
  etc/settings.yaml \
  etc/sample_sheet.csv.example \
  images/Logo_PiGx.png

dist_doc_DATA = 									\
  README.md

EXTRA_DIST = \
  META.in \
  tests/sample_data/Test_vep_sarscov2_parsed.txt \
  tests/sample_data/test_vep_parsing_input.txt \
  tests/sample_sheet.csv \
  tests/settings.yaml \
  $(wildcard databases/sigmut_db/*)

AM_TESTS_ENVIRONMENT = srcdir="$(abs_top_srcdir)" builddir="$(abs_top_builddir)" PIGX_UNINSTALLED=1 PIGX_UGLY=1

TESTS = \
  tests/tests/vep.sh

clean-local:
	-rm -rf .snakemake

CLEANFILES = $(nodist_bin_SCRIPTS) config.json		

integration:
	$(AM_TESTS_ENVIRONMENT) $(abs_top_srcdir)/pigx-sars-cov2-ww --printshellcmds -s tests/settings.yaml tests/sample_sheet.csv

# Build docker image with Guix
docker: dist
	guix pack -C none -e '(load "guix.scm")' --with-source=pigx_sars-cov2-ww-$(VERSION).tar.gz -f docker -S /bin=bin -S /lib=lib -S /share=share glibc-utf8-locales tzdata coreutils bash
