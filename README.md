# CRISPR-DART (Downstream Analysis and Reporting Tool)

crispr-DART is a pipeline to process, analyse, and report about the 
CRISPR-Cas9 induced genome editing outcomes from high-throughput sequencing
of target regions of interest. 

The pipeline allows single/paired-end Illumina reads or long PacBio reads from 
both DNA and RNA samples. 

The pipeline consists of the following steps:
- Quality control (fastqc/multiqc) and improvement (TrimGalore!) of raw reads 
- Mapping the reads to the genome of interest (BBMap)
- Re-alignment of reads with insertions/deletions (GATK)
- Extracting statistics about the detected insertions and deletions
(various R libraries including GenomicAlignments and RSamtools)
- Reporting of the editing outcomes in interactive reports organized into a 
website. (rmarkdown::render_site) 

## Installation

1. Download the source code:

> git clone https://github.com/BIMSBbioinfo/crispr_DART.git

2. Install R/Bioconductor packages

3. Use [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) to install the remaining dependencies

- Create an isolated Conda environment with dependencies
> conda create -n crispr_dart --file requirements.txt

- Activate the environment
> source activate crispr_dart



