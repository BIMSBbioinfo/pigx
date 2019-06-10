<a name="logo"/>
<div align="center">
<img src="images/Logo_PiGx.png" alt="PiGx Logo"  width="30%" height="30%" ></img>
</a>
</div>

**Copyright 2017-2018: Bora Uyar, Jona Ronen, Ricardo Wurmus.**
**This work is distributed under the terms of the GNU General Public License, version 3 or later.  It is free to use for all purposes.**

-----------

## Summary

PiGX RNAseq is an analysis pipeline for preprocessing and reporting for RNA sequencing experiments. It is easy to use and produces high quality reports. The inputs are reads files from the sequencing experiment, and a configuration file which describes the experiment. In addition to quality control of the experiment, the pipeline produces a differential expression report comparing samples in an easily configurable manner.

## What does it do

- Trim reads using trim-galore
- Quality control reads using fastQC and multiQC
- Map reads and quantify read counts per gene using STAR
- Estimate read counts per transcript using SALMON
- Run differential expression analyses using DESeq2

## What does it output

- QC reports
- bam files
- bigwig files
- reads matrix
- DE reports

# Install

You can install this pipeline with all its dependencies using GNU Guix:

    guix package -i pigx-rnaseq

You can also install it from source manually.  You can find the [latest
release](https://github.com/BIMSBbioinfo/pigx_rnaseq/releases/latest)
here.  PiGx uses the GNU build system.  Please make sure that all
required dependencies are installed and then follow these steps after
unpacking the latest release tarball:

```sh
./configure --prefix=/some/where
make install
```

# Dependencies

By default the `configure` script expects tools to be in a directory
listed in the `PATH` environment variable.  If the tools are installed
in a location that is not on the `PATH` you can tell the `configure`
script about them with variables.  Run `./configure --help` for a list
of all variables and options.

You can prepare a suitable environment with Conda or with [GNU
Guix](https://gnu.org/s/guix).  If you do not use one of these package
managers, you will need to ensure that the following software is
installed:

<details>
<summary>Software dependencies</summary>

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
	- rtracklayer
	- SummarizedExperiment
	- gProfileR
- python
	- snakemake
	- pyyaml
- fastqc
- multiqc
- star
- trim-galore
- bedtools
- samtools
- htseq-count

</details>

## Via Conda

Although we highly recommend using Guix to install the software, it is also possible 
to install the dependencies via Conda. 

- Download pigx_rnaseq source code 
    - run: 
    > git clone https://github.com/BIMSBbioinfo/pigx_rnaseq.git
- Download and install Anaconda from https://www.anaconda.com/download
- Locate the 'requirements.txt' file in the source code. 
    - run:
    > conda create --name pigx_rnaseq --file requirements.txt
    - activate the environment:
    > source activate pigx_rnaseq 

## Via Guix

Assuming you have Guix installed, the following command spawns a
sub-shell in which all dependencies are available:

```sh
guix environment -l guix.scm
```


# Getting started

To run PiGx on your experimental data, first enter the necessary parameters in the spreadsheet file (see following section), and then from the terminal type

```sh
$ pigx-rnaseq [options] sample_sheet.csv
```
To see all available options type the `--help` option

```sh
$ pigx-rnaseq --help

usage: pigx-rnaseq [-h] [-v] -s SETTINGS [-c CONFIGFILE] [--target TARGET]
                   [-n] [--graph GRAPH] [--force] [--reason] [--unlock]
                   samplesheet

PiGx RNAseq Pipeline.

PiGx RNAseq is a data processing pipeline for RNAseq read data.

positional arguments:
  samplesheet                             The sample sheet containing sample data in CSV format.

optional arguments:
  -h, --help                              show this help message and exit
  -v, --version                           show program's version number and exit
  -s SETTINGS, --settings SETTINGS        A YAML file for settings that deviate from the defaults.
  -c CONFIGFILE, --configfile CONFIGFILE  The config file used for calling the underlying snakemake process.  By
                                          default the file 'config.json' is dynamically created from the sample
                                          sheet and the settings file.
  --target TARGET                         Stop when the named target is completed instead of running the whole
                                          pipeline.  The default target is "final-report".  Pass "--target=help"
                                          to describe all available targets.
  -n, --dry-run                           Only show what work would be performed.  Do not actually run the
                                          pipeline.
  --graph GRAPH                           Output a graph in Graphviz dot format showing the relations between
                                          rules of this pipeline.  You must specify a graph file name such as
                                          "graph.pdf".
  --force                                 Force the execution of rules, even though the outputs are considered
                                          fresh.
  --reason                                Print the reason why a rule is executed.
  --unlock                                Recover after a snakemake crash.

This pipeline was developed by the Akalin group at MDC in Berlin in 2017-2018.
```


# The input parameters


## Sample Sheet

The sample sheet is a tabular file describing the experiment. The table has the following columns:

| name | reads | reads2 | sample_type | covariate1 |
|------|-------|--------|-------------|------------|

- _name_ is the name for the sample
- _reads1/2_ are the fastq file names of paired end reads
  - the location of these files is specified in `settings.yaml`
  - for single-end data, leave the reads2 column in place, but have it empty
- _sample_type_ can be anything. For instance, a group of biologial replicates.

Additional columns may be included which may be used as covariates in the differential expression analysis (sex, age, different treatments).

## Settings File

The settings file is a _YAML_ file which specifies:

- Locations:
  - The locations of the reads (directory where `fastq` files are)
  - The location of the outputs for the pipeline
  - The location of the `fasta` file with the reference genome (must be prepared by the user)
  - The locations of the transcriptome assembly (for alignment with salmon)
  - The location of a `GTF` file with genome annotations
- Organism (for GO-term analysis using `gProfileR`)
- Differential Expression analyses to be run
  - Which samples to compare (by `sample_type` in the sample sheet)
  - Which covariates to include in the DE analysis (from additional columns in the sample sheet)

In order to get started, enter `pigx-rnaseq --init-settings my_settings.yaml`. This will create a file called `my_settings.yaml` with the default structure. The file will look like this:

```
locations:
  reads-dir: sample_data/reads/
  output-dir: output/
  genome-fasta: sample_data/sample.fasta
  cdna-fasta: sample_data/sample.cdna.fasta
  gtf-file: sample_data/sample.gtf

organism: hsapiens

DEanalyses:
  #names of analyses can be anything but they have to be unique for each combination of case control group comparisons.
  analysis1:
    #if multiple sample names are provided, they must be separated by comma
    case_sample_groups: "HBR"
    control_sample_groups: "UHR"
    covariates: ''

execution:
  submit-to-cluster: no
  jobs: 6
  nice: 19
```

### DEanalysis

The section named `DEanalyses` in the settings file allows the user to specify a number of differential expression analyses to be performed. In the example above, `analysis1` will be the name of the only analysis specified. In that analysis, samples with `sample_type` _HBR_ will be compared with those with `sample_type` _UHR_, with no covariates included in the analysis.

#### Using multiple `sample_type`s as cases or controls

The user may include more than one `sample_type` as the controls or cases for any DE analysis by specifiying a comma-spearated list. For example, the following block specified an analysis which compares samples belonging to _mut1_ and _mut2_ to the _WT_ samples:

```
DEanalyses:
  two_cases_analysis:
    case_sample_groups: "mut1,mut2"
    case_control_groups: "WT"
```

The same may be done to specify several `sample_type`s for the controls.

#### Covariates in DE analysis

Any number of additional columns may be added to the sample sheet and used as covariates in the DE analysis. The following sample sheet includes a column indicating the sex of the sample:

| name       | reads               | reads2              | sample_type | sex |
|------------|---------------------|---------------------|-------------|-----|
| treatment1 | treatment1.r1.fastq | treatment1.r2.fastq | treatment   | m   |
| treatment2 | treatment2.r1.fastq | treatment2.r2.fastq | treatment   | m   |
| control1   | control1.r1.fastq   | control1.r2.fastq   | control     | m   |
| control2   | control2.r1.fastq   | control2.r2.fastq   | control     | f   |

With the use of the following block in the settings file, the sex will be used as a covariate in the differential expression analysis:

```
DEanalyses:
  analysis_with_covariate:
    case_sample_groups: "treatment"
    case_control_groups: "control"
    covariates: "sex"
```

Multiple covariates may be specified by providing a comma-separated list, such as

```
covariates: "sex,age,smoking_history"
```

### Execution

The `execution` section in the settings file allows the user to specify whether the pipeline is to be submitted to a cluster, or run locally, and the degree of parallelism. For a full list of possible parameters, see `etc/settings.yaml`.

# Example

An example can be found in the `tests` directory.  The
`sample_sheet.csv` file here specifies the following sample data:

  - 3 replicates of mRNA from human (UHR)
  - 3 replicates of mRNA from human brain (HBR)
  - only chromosome 22 and ERCC spike-ins

----------------------------------------
2017
