<a name="logo"/>
<div align="center">
<img src="images/Logo_PiGx.png" alt="PiGx Logo"  width="30%" height="30%" ></img>
</a>
</div>

**Copyright 2017-2018: Altuna Akalin.**
**This work is distributed under the terms of the GNU General Public License, version 3 or later.  It is free to use for all purposes.**

-----------

## Summary

PiGX scRNAseq is an analysis pipeline for preprocessing and quality control for single cell RNA sequencing experiments. 
The inputs are read files from the sequencing experiment, and a configuration file which describes the experiment. 
It produces processed files for downstream analysis and interactive quality reports. 
The pipeline is designed to work with UMI based methods. It currently supports all methods which output paired
adapter - read files.


## What does it do

- Quality control reads using fastQC and multiQC
- Automatically determines the appropriate cell number
- Constructs the digital gene expression matrix
- Calculates per sample and per cell statistics
- Prepares a quality control report
- Normalizes data and does dimensionallity reduction


## What does it output

- bam files
- bigwig files
- UMI and read count matrices
- Quality control report
- SingleCellExperiment object with pre-calculated statistics and dimensionallity reductions

# Install

At this time there are no ready-made packages for this pipeline, so
you need to install PiGx from source.

You can find the [latest
release](https://github.com/BIMSBbioinfo/pigx_scrnaseq/releases/latest)
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
Guix](https://gnu.org/s/guix).  

## Via Conda

- Download pigx_scrnaseq source code 
    - run: 
    > git clone https://github.com/BIMSBbioinfo/pigx_scrnaseq.git
- Download and install Anaconda from https://www.anaconda.com/download
- Locate the 'environment.yml' file in the source code. 
    - run:
    > conda env create -f environment.yml #provide path to the environment.yml file
    - activate the environment:
    > source activate pigx_scrnaseq 

## Via Guix

Assuming you have Guix installed, the following command spawns a
sub-shell in which all dependencies are available:

```sh
guix environment -l guix.scm
```

If you do not use one of these package
managers, you will need to ensure that the following software is
installed:

<details>
<summary>Software dependencies</summary>

</details>

# Getting started

To run PiGx on your experimental data, first enter the necessary parameters in the spreadsheet file (see following section), and then from the terminal type.
To run the pipeline, you will also need the appropriate genome sequence in fasta format, and the genome annotation in a 
gtf format.

```sh
$ pigx-rscnaseq [options] sample_sheet.csv -s settings.yaml
```

To see all available options type the `--help` option

```sh
$ pigx-scrnaseq --help

usage: pigx-scrnaseq [-h] [-v] -s SETTINGS [-c CONFIGFILE] [--target TARGET]
                   [-n] [--graph GRAPH] [--force] [--reason] [--unlock]
                   samplesheet

PiGx scRNAseq Pipeline.

PiGx scRNAseq is a data processing pipeline for single cell RNAseq read data.

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

| name | reads1 | reads2 | library | covariate1 | covariate2 |
|------|--------|--------|---------|------------|------------|

- _name_ - name for the sample, which will be used to label the sample in all downstream analysis
- _reads1 - fastq file containing the **adapter sequences**
- _reads2 - fastq file containing the **sequenced reads**
  - location of these files is specified in `settings.yaml`
- _library - sequencing platform on which the experiment was performed (i.e. dropseq)
- _covariates - variables which describe the samples

Additional columns may be included which may be used as covariates in the differential expression analysis (sex, age, different treatments).

## Settings File

The settings file is a _YAML_ file which specifies:

- Locations:
  - The locations of the reads (directory where `fastq` files are located)
  - The location of the output directory
  - The location of the `fasta` file with the reference genome (must be prepared by the user)
  - The location of a `GTF` file with genome annotations
- Organism name

In order to get started, enter `pigx-scrnaseq --init-settings my_settings.yaml`. This will create a file called `my_settings.yaml` with the default structure. The file will look like this:

```
locations:
  reads-dir: sample_data/reads/
  output-dir: output/
  genome-fasta: sample_data/sample.fasta
  cdna-fasta: sample_data/sample.cdna.fasta
  gtf-file: sample_data/sample.gtf

organism: hsapiens

execution:
  submit-to-cluster: no
  jobs: 6
  nice: 19
```

### Mixed species experiments

**TODO**

# Output file description

# Detailed pipeline description

# Downstream analysis

# Using iSEE for interactive exploration of the single cell experiment object

### Cluster Execution

The `execution` section in the settings file allows the user to specify whether the pipeline is to be submitted to a cluster, or run locally, and the degree of parallelism. For a full list of possible parameters, see `etc/settings.yaml`.

# Example

An example can be found in the `tests` directory.  The
`sample_sheet.csv` file here specifies the following sample data:

 
----------------------------------------
2018
