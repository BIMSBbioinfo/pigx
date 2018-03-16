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
  scripts/counts_matrix_from_SALMON.R		\
  scripts/validate_input.py

dist_pkgdata_DATA =									\
  etc/sample_sheet.csv.example								\
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
  tests/sample_data/sample.cdna.fasta				\
  tests/sample_data/preprocessed/trimmed_reads/HBR_Rep1.read2.fastq.gz_trimming_report.txt \
  tests/sample_data/preprocessed/trimmed_reads/UHR_Rep3.read2.fastq.gz_trimming_report.txt \
  tests/sample_data/preprocessed/trimmed_reads/HBR_Rep2.read1.fastq.gz_trimming_report.txt \
  tests/sample_data/preprocessed/trimmed_reads/UHR_Rep3_R2.fastq.gz \
  tests/sample_data/preprocessed/trimmed_reads/UHR_Rep3.read1.fastq.gz_trimming_report.txt \
  tests/sample_data/preprocessed/trimmed_reads/HBR_Rep3_R2.fastq.gz \
  tests/sample_data/preprocessed/trimmed_reads/UHR_Rep1.read1.fastq.gz_trimming_report.txt \
  tests/sample_data/preprocessed/trimmed_reads/UHR_Rep1_R1.fastq.gz \
  tests/sample_data/preprocessed/trimmed_reads/UHR_Rep2.read2.fastq.gz_trimming_report.txt \
  tests/sample_data/preprocessed/trimmed_reads/HBR_Rep2_R2.fastq.gz \
  tests/sample_data/preprocessed/trimmed_reads/UHR_Rep2_R1.fastq.gz \
  tests/sample_data/preprocessed/trimmed_reads/HBR_Rep3.read2.fastq.gz_trimming_report.txt \
  tests/sample_data/preprocessed/trimmed_reads/UHR_Rep2_R2.fastq.gz \
  tests/sample_data/preprocessed/trimmed_reads/HBR_Rep3.read1.fastq.gz_trimming_report.txt \
  tests/sample_data/preprocessed/trimmed_reads/HBR_Rep3_R1.fastq.gz \
  tests/sample_data/preprocessed/trimmed_reads/UHR_Rep1_R2.fastq.gz \
  tests/sample_data/preprocessed/trimmed_reads/UHR_Rep3_R1.fastq.gz \
  tests/sample_data/preprocessed/trimmed_reads/HBR_Rep1_R2.fastq.gz \
  tests/sample_data/preprocessed/trimmed_reads/HBR_Rep1_R1.fastq.gz \
  tests/sample_data/preprocessed/trimmed_reads/HBR_Rep2.read2.fastq.gz_trimming_report.txt \
  tests/sample_data/preprocessed/trimmed_reads/HBR_Rep2_R1.fastq.gz \
  tests/sample_data/preprocessed/trimmed_reads/UHR_Rep1.read2.fastq.gz_trimming_report.txt \
  tests/sample_data/preprocessed/trimmed_reads/HBR_Rep1.read1.fastq.gz_trimming_report.txt \
  tests/sample_data/preprocessed/trimmed_reads/UHR_Rep2.read1.fastq.gz_trimming_report.txt \
  tests/sample_data/preprocessed/colData.tsv \
  tests/sample_data/preprocessed/fastqc/HBR_Rep3_Aligned.sortedByCoord.out_fastqc.zip \
  tests/sample_data/preprocessed/fastqc/UHR_Rep1_Aligned.sortedByCoord.out_fastqc.zip \
  tests/sample_data/preprocessed/fastqc/UHR_Rep2_Aligned.sortedByCoord.out_fastqc.zip \
  tests/sample_data/preprocessed/fastqc/UHR_Rep3_Aligned.sortedByCoord.out_fastqc.zip \
  tests/sample_data/preprocessed/fastqc/HBR_Rep1_Aligned.sortedByCoord.out_fastqc.zip \
  tests/sample_data/preprocessed/fastqc/HBR_Rep2_Aligned.sortedByCoord.out_fastqc.zip \
  tests/sample_data/preprocessed/preprocessed_data/counts_from_STAR.tsv \
  tests/sample_data/preprocessed/mapped_reads/UHR_Rep3_Log.out \
  tests/sample_data/preprocessed/mapped_reads/HBR_Rep2_Aligned.sortedByCoord.out.bam \
  tests/sample_data/preprocessed/mapped_reads/UHR_Rep1_Log.progress.out \
  tests/sample_data/preprocessed/mapped_reads/HBR_Rep2_Aligned.toTranscriptome.out.bam \
  tests/sample_data/preprocessed/mapped_reads/HBR_Rep1_Log.progress.out \
  tests/sample_data/preprocessed/mapped_reads/HBR_Rep2_Log.final.out \
  tests/sample_data/preprocessed/mapped_reads/HBR_Rep1_ReadsPerGene.out.tab \
  tests/sample_data/preprocessed/mapped_reads/HBR_Rep3_Aligned.sortedByCoord.out.bam \
  tests/sample_data/preprocessed/mapped_reads/UHR_Rep2_SJ.out.tab \
  tests/sample_data/preprocessed/mapped_reads/UHR_Rep3_ReadsPerGene.out.tab \
  tests/sample_data/preprocessed/mapped_reads/HBR_Rep1_Log.final.out \
  tests/sample_data/preprocessed/mapped_reads/UHR_Rep3_Log.progress.out \
  tests/sample_data/preprocessed/mapped_reads/HBR_Rep1_Aligned.sortedByCoord.out.bam \
  tests/sample_data/preprocessed/mapped_reads/UHR_Rep2_Log.final.out \
  tests/sample_data/preprocessed/mapped_reads/HBR_Rep1_Log.out \
  tests/sample_data/preprocessed/mapped_reads/UHR_Rep1_Log.out \
  tests/sample_data/preprocessed/mapped_reads/HBR_Rep3_Log.final.out \
  tests/sample_data/preprocessed/mapped_reads/UHR_Rep2_Log.progress.out \
  tests/sample_data/preprocessed/mapped_reads/HBR_Rep2_SJ.out.tab \
  tests/sample_data/preprocessed/mapped_reads/UHR_Rep2_Log.out \
  tests/sample_data/preprocessed/mapped_reads/UHR_Rep1_ReadsPerGene.out.tab \
  tests/sample_data/preprocessed/mapped_reads/HBR_Rep2_Log.progress.out \
  tests/sample_data/preprocessed/mapped_reads/UHR_Rep3_Aligned.sortedByCoord.out.bam \
  tests/sample_data/preprocessed/mapped_reads/HBR_Rep1_SJ.out.tab \
  tests/sample_data/preprocessed/mapped_reads/HBR_Rep2_ReadsPerGene.out.tab \
  tests/sample_data/preprocessed/mapped_reads/UHR_Rep3_Log.final.out \
  tests/sample_data/preprocessed/mapped_reads/HBR_Rep3_Log.out \
  tests/sample_data/preprocessed/mapped_reads/UHR_Rep3_Aligned.toTranscriptome.out.bam \
  tests/sample_data/preprocessed/mapped_reads/UHR_Rep1_Log.final.out \
  tests/sample_data/preprocessed/mapped_reads/UHR_Rep1_SJ.out.tab \
  tests/sample_data/preprocessed/mapped_reads/UHR_Rep1_Aligned.sortedByCoord.out.bam \
  tests/sample_data/preprocessed/mapped_reads/UHR_Rep2_ReadsPerGene.out.tab \
  tests/sample_data/preprocessed/mapped_reads/HBR_Rep3_ReadsPerGene.out.tab \
  tests/sample_data/preprocessed/mapped_reads/HBR_Rep2_Log.out \
  tests/sample_data/preprocessed/mapped_reads/UHR_Rep3_SJ.out.tab \
  tests/sample_data/preprocessed/mapped_reads/HBR_Rep3_Log.progress.out \
  tests/sample_data/preprocessed/mapped_reads/UHR_Rep2_Aligned.toTranscriptome.out.bam \
  tests/sample_data/preprocessed/mapped_reads/UHR_Rep2_Aligned.sortedByCoord.out.bam \
  tests/sample_data/preprocessed/mapped_reads/HBR_Rep1_Aligned.toTranscriptome.out.bam \
  tests/sample_data/preprocessed/mapped_reads/HBR_Rep3_Aligned.toTranscriptome.out.bam \
  tests/sample_data/preprocessed/mapped_reads/HBR_Rep3_SJ.out.tab \
  tests/sample_data/preprocessed/mapped_reads/UHR_Rep1_Aligned.toTranscriptome.out.bam \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep3/quant.genes.sf \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep3/quant.sf \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep3/logs/salmon_quant.log \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep3/aux_info/obs_gc.gz \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep3/aux_info/exp_gc.gz \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep3/aux_info/obs5_seq.gz \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep3/aux_info/observed_bias.gz \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep3/aux_info/observed_bias_3p.gz \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep3/aux_info/exp5_seq.gz \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep3/aux_info/ambig_info.tsv \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep3/aux_info/meta_info.json \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep3/aux_info/expected_bias.gz \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep3/aux_info/fld.gz \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep3/aux_info/exp3_seq.gz \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep3/aux_info/obs3_seq.gz \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep3/cmd_info.json \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep3/libParams/flenDist.txt \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep3/lib_format_counts.json \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep2/quant.genes.sf \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep2/quant.sf \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep2/logs/salmon_quant.log \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep2/aux_info/obs_gc.gz \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep2/aux_info/exp_gc.gz \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep2/aux_info/obs5_seq.gz \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep2/aux_info/observed_bias.gz \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep2/aux_info/observed_bias_3p.gz \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep2/aux_info/exp5_seq.gz \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep2/aux_info/ambig_info.tsv \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep2/aux_info/meta_info.json \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep2/aux_info/expected_bias.gz \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep2/aux_info/fld.gz \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep2/aux_info/exp3_seq.gz \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep2/aux_info/obs3_seq.gz \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep2/cmd_info.json \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep2/libParams/flenDist.txt \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep2/lib_format_counts.json \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep2/quant.genes.sf \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep2/quant.sf \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep2/logs/salmon_quant.log \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep2/aux_info/obs_gc.gz \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep2/aux_info/exp_gc.gz \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep2/aux_info/obs5_seq.gz \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep2/aux_info/observed_bias.gz \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep2/aux_info/observed_bias_3p.gz \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep2/aux_info/exp5_seq.gz \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep2/aux_info/ambig_info.tsv \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep2/aux_info/meta_info.json \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep2/aux_info/expected_bias.gz \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep2/aux_info/fld.gz \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep2/aux_info/exp3_seq.gz \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep2/aux_info/obs3_seq.gz \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep2/cmd_info.json \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep2/libParams/flenDist.txt \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep2/lib_format_counts.json \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep1/quant.genes.sf \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep1/quant.sf \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep1/logs/salmon_quant.log \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep1/aux_info/obs_gc.gz \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep1/aux_info/exp_gc.gz \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep1/aux_info/obs5_seq.gz \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep1/aux_info/observed_bias.gz \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep1/aux_info/observed_bias_3p.gz \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep1/aux_info/exp5_seq.gz \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep1/aux_info/ambig_info.tsv \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep1/aux_info/meta_info.json \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep1/aux_info/expected_bias.gz \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep1/aux_info/fld.gz \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep1/aux_info/exp3_seq.gz \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep1/aux_info/obs3_seq.gz \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep1/cmd_info.json \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep1/libParams/flenDist.txt \
  tests/sample_data/preprocessed/salmon_output/HBR_Rep1/lib_format_counts.json \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep3/quant.genes.sf \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep3/quant.sf \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep3/logs/salmon_quant.log \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep3/aux_info/obs_gc.gz \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep3/aux_info/exp_gc.gz \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep3/aux_info/obs5_seq.gz \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep3/aux_info/observed_bias.gz \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep3/aux_info/observed_bias_3p.gz \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep3/aux_info/exp5_seq.gz \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep3/aux_info/ambig_info.tsv \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep3/aux_info/meta_info.json \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep3/aux_info/expected_bias.gz \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep3/aux_info/fld.gz \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep3/aux_info/exp3_seq.gz \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep3/aux_info/obs3_seq.gz \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep3/cmd_info.json \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep3/libParams/flenDist.txt \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep3/lib_format_counts.json \
  tests/sample_data/preprocessed/salmon_output/counts_from_SALMON.tsv \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep1/quant.genes.sf \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep1/quant.sf \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep1/logs/salmon_quant.log \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep1/aux_info/obs_gc.gz \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep1/aux_info/exp_gc.gz \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep1/aux_info/obs5_seq.gz \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep1/aux_info/observed_bias.gz \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep1/aux_info/observed_bias_3p.gz \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep1/aux_info/exp5_seq.gz \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep1/aux_info/ambig_info.tsv \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep1/aux_info/meta_info.json \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep1/aux_info/expected_bias.gz \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep1/aux_info/fld.gz \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep1/aux_info/exp3_seq.gz \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep1/aux_info/obs3_seq.gz \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep1/cmd_info.json \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep1/libParams/flenDist.txt \
  tests/sample_data/preprocessed/salmon_output/UHR_Rep1/lib_format_counts.json \
  tests/test_trim_galore/input_files/sample.read1.fastq.gz \
  tests/test_trim_galore/input_files/sample.read2.fastq.gz \
  tests/settings.yaml \
  tests/settings_no_de.yaml

AM_TESTS_ENVIRONMENT = srcdir="$(abs_top_srcdir)" builddir="$(abs_top_builddir)" PIGX_UNINSTALLED=1 PIGX_UGLY=1

TESTS = \
  tests/test_genome_coverage/test.sh					\
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
	guix pack -C none -e '(load "guix.scm")' --with-source=pigx_rnaseq-$(VERSION).tar.gz -f docker -S /bin=bin -S /lib=lib -S /share=share glibc-utf8-locales tzdata coreutils bash
