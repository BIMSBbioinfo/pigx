# pigx pipeline for RNAseq



<a name="logo"/>
<div align="center">
<img src="images/Logo_PiGx.png" alt="PiGx Logo"  width="30%" height="30%" ></img>
</a>
</div>

**Copyright 2017: Bora Uyar, Jona Ronen, Ricardo Wurmus.**
**This work is distributed under the terms of the GNU General Public License, version 3 or later.  It is free to use for all purposes.**


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

# Install

At this time there are no ready-made packages for this pipeline, so
you need to install PiGx from source.

You can find the [latest
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
- bamCoverage
- samtools
- htseq-count


## Via Conda

- Download pigx_rnaseq source code 
    - run: 
    > git clone https://github.com/BIMSBbioinfo/pigx_rnaseq.git
- Download and install Anaconda from https://www.anaconda.com/download
- Locate the 'environment.yml' file in the source code. 
    - run:
    > conda env create -f environment.yml #provide path to the environment.yml file
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

usage: pigx-rnaseq [-h] [-v] [-s SETTINGS] [-c CONFIGFILE] [--snakeparams SNAKEPARAMS]
                   samplesheet

PiGx RNAseq Pipeline.

PiGx is a data processing pipeline for RNA seq read data.

positional arguments:
  samplesheet                             The sample sheet containing sample data in CSV format.

optional arguments:
  -h, --help                              show this help message and exit
  -v, --version                           show program's version number and exit
  -s SETTINGS, --settings SETTINGS        A YAML file for settings that deviate from the defaults.
  -c CONFIGFILE, --configfile CONFIGFILE  The config file used for calling the underlying snakemake process.  By
                                          default the file 'config.json' is dynamically created from the sample
                                          sheet and the settings file.
  --snakeparams SNAKEPARAMS               Additional parameters to be passed down to snakemake, e.g.
                                              --dryrun    do not execute anything
                                              --forceall  re-run the whole pipeline
```


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
  - for single-end data, leave the reads2 column in place, but have it empty
- sample_type can be anything
- comparison_factor should say which samples to compare together for DE analysis
  - so samples that are to be compared together should have the same factor

## Example

An example can be found in the `tests` directory.  The
`sample_sheet.csv` file here specifies the following sample data:

  - 3 replicates of mRNA from human (UHR)
  - 3 replicates of mRNA from human brain (HBR)
  - only chromosome 22 and ERCC spike-ins

## To do

- fix gene sets for DE report

----------------------------------------
2017
