

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


## Installation of dependecies via Guix


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
$ pigx_bsseq [options] tablesheet
```
To see all available options type the `--help` option
```sh
$ pigx_bsseq --help

usage: pigx_bsseq [-h] [-v] [-p PROGRAMS] [-c CONFIGFILE] [-s SNAKEPARAMS]
               tablesheet

PiGx BSseq Pipeline.

PiGx is a data processing pipeline for raw fastq read data of
bisulfite experiments.  It produces methylation and coverage
information and can be used to produce information on differential
methylation and segmentation.

positional arguments:
  tablesheet                                 The tablesheet containing the basic configuration information for
                                             running the pipeline.

optional arguments:
  -h, --help                                 show this help message and exit
  -v, --version                              show program's version number and exit
  -p PROGRAMS, --programs PROGRAMS           A JSON file containing the absolute paths of the required tools.
  -c CONFIGFILE, --configfile CONFIGFILE     The config file used for calling the underlying snakemake process.  By
                                             default the file 'config.json' is dynamically created from tablesheet
                                             and programs file.
  -s SNAKEPARAMS, --snakeparams SNAKEPARAMS  Additional parameters to be passed down to snakemake, e.g.
                                                 --dryrun    do not execute anything
                                                 --forceall  re-run the whole pipeline
```

# Input parameters

The input parameters specifying the desired behaviour of PiGx should
be entered into the tablesheet file.  When PiGx is run, the data from
this file will be used to automatically generate a configuration file.

Here is an example tablesheet:
```
[ GENERAL PARAMETERS ]
PATHIN="in/"
PATHOUT="out/"
GENOMEPATH="genome/"
GENOME_VERSION="hg19"
bismark_args=" -N 0 -L 20 "
fastqc_args=""
trim_galore_args=""
bam_methCall_args_mincov="0"
bam_methCall_args_minqual="10"
NICE="19"
numjobs="6"
cluster_run="FALSE"
contact_email="NONE"
bismark_cores="3"
bismark_MEM="19G"
MEM_default="8G"
qname="all"
h_stack="128m"
diffmeth_cores="20"


[ SAMPLES ]
Read1,Read2,SampleID,ReadType,Treatment
PE_1.fq.gz,PE_2.fq.gz,PEsample,WGBS,0
SE_techrep1.fq.gz,,SEsample,WGBS,1
SE_techrep2.fq.gz,,SEsample_v2,WGBS,2

[ DIFFERENTIAL METHYLATION ]
0, 1
```

The tablesheet contains 3 paragraphs: 
- general parameters, 
- a table with sample specific information containing the names of fastq files, unique sample ids, the type of bisulfite sequencing experiment (could be RRBS or WGBS,only WGBS is available right now) and treatment group for differential methylation detection
- treatment groups considered for differential methylation detection

## Details about General Parameters

General parameters have to contain variables:

<details>
  <summary>Click to expand explanations</summary>

| Variable name | description |
| ------------- |:-----------:|
| PATHIN        | string: location of the experimental\nall input data files (.fastq[.gz\|.bz2])   |
| PATHOUT       | string: ultimate location of the output data and report files   |
| GENOMEPATH    | string: location of the reference genome data for alignment   |
| GENOME_VERSION| string: an UCSC assembly release name e.g. "hg19"
| bismark_args  | string: optional arguments supplied to bismark during alignment. See the [Bismark User Guide], e.g. " -N 0 -L 20 "
| fastqc_args  | string: optional arguments supplied to FastQC during alignment. See the [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), e.g. ""
| trim_galore_args | string: optional arguments supplied to Trim Galore! during alignment. See the [Trim Galore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) e.g. "" 
| bam_methCall_args_mincov | string: minimum read coverage to be included in the methylKit objects. defaults to 10. Any methylated base/region in the text files below the mincov value will be ignored.
| bam_methCall_args_minqual | string: minimum phred quality score to call a methylation status for a base, e.g. "10"
| cluster_run | string: a boolean whether the pipeline should be run on cluster, e.g. "FALSE"
| numjobs | string: number of jobs sent to cluster, e.g. "6"
| contact_email | string: email address to which information about cluster job is sent
| bismark_cores | string: number of cores used by bismark, e.g. "3"
| bismark_MEM | string: amount of memory used by bismark, e.g. "19G"
| MEM_default | string: amount of memory used for all jobs besides bismark, e.g. "8G"
| qname | string: queue name (used for cluster jobs), e.g. "all"
| h_stack | string: stack size limit (used for cluster jobs), e.g. "128m"
| diffmeth_cores | integer: denoting how many cores should be used for parallel differential methylation calculations
| NICE          | integer: from -20 to 19; higher values make the program execution less demanding on computational resources 

</details>
 
 
Make sure that all input files (paired or single end) are present in the folder
indicated by `PATHIN`. All output produced by the pipeline will written to the folder indicated by `PATHOUT`,
with subdirectories corresponding to the various stages of the process.
The directory pointed to by `GENOMEPATH` has to contain the reference genome being mapped to.
