

<a name="logo"/>
<div align="center">
<img src="images/Logo_PiGx.png" alt="PiGx Logo"  width="30%" height="30%" ></img>
</a>
</div>

**Copyright 2017: Alexander Gosdschan, Katarzyna Wreczycka, Bren Osberg, Ricardo Wurmus.**
**This work is distributed under the terms of the GNU General Public License, version 3 or later.  It is free to use for all purposes.**

-----
# Summary

PiGx is a data processing pipeline for raw fastq read data of bisulfite experiments; it produces reports on aggregate methylation and  coverage and can be used to produce information on differential methylation and segmentation. It was first developed by the Akalin group at MDC in Berlin in 2017.

The figure below provides a sketch of the process.
![](images/pipelineIO_BSseq.png )


# Install

PiGx uses the GNU build system.  If you want to install PiGx from
source (here you can find the [latest release](https://github.com/BIMSBbioinfo/pigx_bsseq/releases/latest)), please make sure that all required dependencies are installed and 
then follow these steps after unpacking the latest release tarball:

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

<details>
  <summary>The following tools must be available:</summary>

 - fastqc
 - trim_galore
 - cutadapt
 - bismark_genome_preparation
 - deduplicate_bismark
 - bismark
 - bowtie2
 - samtools [>=1.3]
 - snakemake
 - Python [>=3.5]
 - PyYAML
 - [pandoc](http://pandoc.org/)
 - [pandoc-citeproc](http://pandoc.org/)
 - R
 - [methylKit](https://github.com/al2na/methylKit) [>=1.3.1]
 - [genomation](http://bioinformatics.mdc-berlin.de/genomation/)
 - [GenomeInfoDb](https://www.bioconductor.org/packages/release/bioc/html/GenomeInfoDb.html)
 - [DT](https://rstudio.github.io/DT/) 
 - [annotationhub](https://www.bioconductor.org/packages/release/bioc/html/AnnotationHub.html)
 - [rtracklayer](http://bioconductor.org/packages/release/bioc/html/rtracklayer.html)
 - [rmarkdown](http://rmarkdown.rstudio.com/) [>=1.5]
 - [bookdown](https://github.com/rstudio/bookdown/)

All of these dependencies must be present in the environment at
configuration time.
</details>


## Installation of dependencies via Guix


You can install PiGx through Guix (TODO: add details here after release).

Run the `configure` script to probe your environment for tools needed
by the pipeline.  If you cannot be bothered to install all packages
manually, we recommend using [GNU Guix](https://gnu.org/s/guix).  The
following command spawns a sub-shell in which all dependencies are
available:

```sh
guix environment -l guix.scm
```

# Getting started

To run PiGx on your experimental data, first enter the necessary parameters in the spreadsheet file (see following section), and then from the terminal type

```sh
$ pigx-bsseq [options] samplesheet.csv
```
To see all available options type the `--help` option
```sh
$ pigx-bsseq --help

usage: pigx-bsseq [-h] [-v] [-s SETTINGS] [-c CONFIGFILE] [--snakeparams SNAKEPARAMS]
                  samplesheet

PiGx BSseq Pipeline.

PiGx is a data processing pipeline for raw fastq read data of
bisulfite experiments.  It produces methylation and coverage
information and can be used to produce information on differential
methylation and segmentation.

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

# Input parameters

The pipeline expects two kinds of input: a sample sheet in CSV format
and a settings file specifying the desired behaviour of PiGx.  PiGx
will automatically generate a JSON configuration file from these
inputs.

The sample sheet is a table with sample-specific information
containing the names of fastq files, unique sample ids, the type of
bisulfite sequencing experiment (could be RRBS or WGBS,only WGBS is
available right now) and treatment group for differential methylation
detection.

Here is an example sample sheet:

```
Read1,Read2,SampleID,ReadType,Treatment
PE_1.fq.gz,PE_2.fq.gz,PEsample,WGBS,0
SE_techrep1.fq.gz,,SEsample,WGBS,1
SE_techrep2.fq.gz,,SEsample_v2,WGBS,2
```

The default settings file can be found at `etc/settings.yaml`.
Settings that are not specified by the user are taken from the default
settings file.  A user's settings file to override some defaults might
look something like this:

```
general:
  methylation-calling:
    minimum-coverage: 1
    minimum-quality: 3
  differential-methylation:
    cores: 2
    treatment-groups:
      - ['A', 'B']

execution:
  submit-to-cluster: yes
  jobs: 6
  cluster:
    memory: 8G
    stack: 128M
    queue: all
    contact-email: foo@example.com

tools:
  bismark:
    cores: 3
```


## Available settings

PiGx recognizes four sections in the settings file:

- `locations` for input, output, and genome directories
- `general` for general settings
- `execution` for settings affecting the pipeline execution
- `tools` for tool-specific paths and arguments.

### Locations

| Variable name | description |
| ------------- |:-----------:|
| input-dir     | string: location of the experimental\nall input data files (.fastq[.gz\|.bz2]) |
| output-dir    | string: ultimate location of the output data and report files   |
| genome-dir    | string: location of the reference genome data for alignment   |

Make sure that all input files (paired or single end) are present in
the folder indicated by `input-dir`. All output produced by the
pipeline will written to the folder indicated by `output-dir`, with
subdirectories corresponding to the various stages of the process.
The directory pointed to by `genome-dir` has to contain the reference
genome being mapped to.

### General

```
general:
  genome-version: hg19
  methylation-calling:
    minimum-coverage: 0
    minimum-quality: 10
  differential-methylation:
    cores: 20
    treatment-groups:
      - ['0', '1']
```

| Variable name | description |
| ------------- |:-----------:|
| genome-version | string: an UCSC assembly release name e.g. "hg19"
| methylation-calling:minimum-coverage | integer: minimum read coverage to be included in the methylKit objects. Defaults to 10. Any methylated base/region in the text files below the mincov value will be ignored.
| methylation-calling:minimum-quality | integer: minimum phred quality score to call a methylation status for a base.  Defaults to 10.
| differential-methylation:cores | integer: denotes how many cores should be used for parallel differential methylation calculations
| differential-methylation:treatment-groups | array of strings

### Execution

```
execution:
  submit-to-cluster: no
  jobs: 6
  nice: 19
  cluster:
    memory: 8G
    stack: 128M
    queue: all
    contact-email: none
```

| Variable name         | description |
| --------------------- |:-----------:|
| submit-to-cluster     | string: whether the pipeline should run locally ("no") or on a cluster
| jobs                  | string: number of jobs sent to cluster, e.g. "6"
| nice                  | integer: from -20 to 19; higher values make the program execution less demanding on computational resources
| cluster:memory        | string: amount of memory used for all jobs besides bismark, e.g. "8G"
| cluster:stack         | string: stack size limit (used for cluster jobs), e.g. "128m"
| cluster:queue         | string: queue name (used for cluster jobs), e.g. "all"
| cluster:contact-email | string: email address to which information about cluster job is sent

### Tools

The values for the `executable` field for each tool are determined at
configure time and usually won't have to be changed unless you want to
experiment with a custom variant of a particular tool.

The `args` field for each tool accepts a string for additional
arguments to be passed to the specified tool.

The `bismark` tool supports additional settings, such as `cores` (the
number of cores used by bismark) and `memory` for the amount of RAM
that bismark may use.
 
