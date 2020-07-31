<p align="center">
  <img alt="logo" src="img/logo.png" width="35%" height="35%">
</p>

# crispr-DART (Downstream Analysis and Reporting Tool)

crispr-DART is a pipeline to process, analyse, and report about the 
CRISPR-Cas9 induced genome editing outcomes from high-throughput sequencing
of target regions of interest. 

crispr-DART has been developed as part of the study "Parallel genetics of regulatory sequences in vivo"
by Froehlich & Uyar et al, 2020 ([pre-print available on Biorxiv](https://www.biorxiv.org/content/10.1101/2020.07.28.224998v2)). 

## Pipeline scheme

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

![pipeline](./img/pipeline_scheme.jpg)


## Example HTML report output

The HTML reports produced by the pipeline are automatically organised as a website. 
Example report website can be browsed [here](https://bimsbstatic.mdc-berlin.de/akalin/buyar/froehlich_uyar_et_al_2020/reports/index.html): 

## Example screenshots from the reports

You can find below some example screenshots from the HTML reports:

![screenshots](./img/reports_screenshots.jpg)


# Installation

1. Download the source code:

```
> git clone https://github.com/BIMSBbioinfo/crispr_DART.git
```

2. Install R/Bioconductor packages

```
> if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

> BiocManager::install(c('data.table', 'yaml', 'ggplot2', 'knitr', 'ggrepel', 'pbapply', 'DT', 
'Biostrings', 'GenomicAlignments', 'rtracklayer', 'GenomicRanges', 'Rsamtools', 'reshape2', 'GenomeInfoDb',
'fastseg', 'gtools', 'IRanges', 'rmarkdown'))
```

3. Use [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) to install the remaining dependencies

- Create an isolated Conda environment with dependencies

```
> conda create -n crispr_dart --file requirements.txt
```

In addition, for the indel re-alignment step, GenomeAnalysisToolKit.jar file for the GATK version 3.8.0 is required. 
The GenomeAnalysisToolKit.jar file needs to be downloaded from 
https://software.broadinstitute.org/gatk/download/archive and stored somewhere that is accessible to the pipeline (e.g. ~/tools/)

- Activate the environment

```
> source activate crispr_dart
```

4. Test the installation

The pipeline can be simply tested by running the bash script `test.sh`. 

The test script uses the necessary input files available in the `sample_data` folder 
and runs the pipeline. If this test runs to completion, you should be ready to analyse your own
data. 


```
> bash ./test.sh
```

# How to run the pipeline 

## Preparing the input files

The pipeline currently requires four different input files. 
1. A sample sheet file, which describes the samples, associated fastq files, the sets of sgRNAs used in the sample and the list of regions of interest. 

Please see the example sample sheet file under `sample_data/sample_sheet.csv`. 

2. A BED file containing the genomic coordinates of all the sgRNAs used in this project. 

Please see the example BED file for sgRNA target sites under `sample_data/cut_sites.bed`

3. A comparisons table, which is used for comparing pairs of samples in terms of genome editing outcomes. 

Please see the example table under `sample_data/comparisons.tsv`

4. A settings file, which combines all the information from the other input files and additional configurations for resource requirements of tools. 

Please see the example file under `sample_data/settings.yaml`

The `sample_data/fasta` folder contains fasta format sequence files that are used as the target genome sequence. 
The `sample_data/reads` folder contains sample read files (fastq.gz files from Illumina and PacBio sequenced samples). 

## Running the pipeline

Once the `settings.yaml` file is configured with paths to all the other required files, the pipeline can simply be run using the bash script `run.sh` requesting 2 cpus. 

```
> bash run.sh */path/to/settings.yaml* 2  
```

If you would like to do a dry-run, meaning that the list of jobs are created but not executed, you can do 


```
> bash run.sh */path/to/settings.yaml* 2 --dry
```

Any additional arguments to `run.sh` after the argument for the number of cpus are passed as arguments to `snakemake`. 

# How to cite

See the pre-print on [Biorxiv](https://www.biorxiv.org/content/10.1101/2020.07.28.224998v2)

# Credits

The software has been developed by [Bora Uyar](https://scholar.google.com/citations?user=YEZr1LUAAAAJ&hl=en) from the [Akalin Lab](https://bioinformatics.mdc-berlin.de) with significant conceptual contributions by [Jonathan Froehlich](https://scholar.google.com/citations?user=aXWiWfcAAAAJ&hl=en) from the [N.Rajewsky Lab](https://www.mdc-berlin.de/n-rajewsky) at the Berlin Institute of Medical Systems Biology of the Max-Delbruck-Center for Molecular Medicine. 

