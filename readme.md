# pigx pipeline for RNAseq

## What does it do

- trim reads
- quality control reads
- map reads using STAR
- count genes using HTseq

## What does it output

- qc report
- bam files
- bigwig files
- reads matrix
- DE report

## How to configure

### `settings.yaml`

- This file points to the folder containing reads, genome index, and annotations
- it also points to where executables of different tools lie

### `sample_sheet.csv`

- has the following columns:

| name | reads | reads2 | sample_type | comparison_factor |
|------|-------|--------|-------------|-------------------|

- name is the name for the sample
- reads1/2 are the fastq file names of paired end reads
  - the location of these files is specified in `settings.yaml`
- sample_type can be anything
- comparison_factor should say which samples to compare together for DE analysis
  - so samples that are to be compared together should have the same factor

## How to run

- once `settings.yaml` and `sample_sheet.csv` have been modified,

	snakemake -s pigx_rnaseq.py

## Example

- the `sample_sheet.csv` file here points to some data on `/data/akalin/`
  - 3 replicates of mRNA from human (UHR)
  - 3 replicates of mRNA from human brain (HBR)
  - only chromosome 22 and ERCC spike-ins

## To do

- configure script for `settings.yaml`
- have pipeline make reference genome
- fix gene sets for DE report
- add support for single end reads

## Software requirements

- R
	- ggplot2
	- ggrepel
	- DESeq2
	- DT
	- pheatmap
	- dendsort
	- corrplot
	- reshape2
	- plotly
	- scales
	- crosstalk
	- gage
- python
	- snakemake
	- pyyaml
	- pandas
- fastqc
- multiqc
- star
- trim-galore
- bamCoverage
- samtools
- htseq-count

----------------------------------------
2017