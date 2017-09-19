

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

You can install PiGx through Guix (TODO: add details here after release).

PiGx uses the GNU build system.  If you want to install PiGx from
source, follow these steps after unpacking the latest release tarball:

```sh
./configure --prefix=/some/where
make install
```

If you are a developer or want to bootstrap the tarball yourself:

```sh
./bootstrap.sh
./configure
make distcheck
```

To run the pipeline without installing it set the environment variable
`PIGX_BSSEQ_UNINSTALLED` before running the pipeline script.


# Dependencies

Run the `configure` script to probe your environment for tools needed
by the pipeline.  If you cannot be bothered to install all packages
manually, we recommend using [GNU Guix](https://gnu.org/s/guix).  The
following command spawns a sub-shell in which all dependencies are
available:

```sh
guix environment -l guix.scm
```

By default the `configure` script expects tools to be in a directory
listed in the `PATH` environment variable.  If the tools are installed
in a location that is not on the `PATH` you can tell the `configure`
script about them with variables.  Run `./configure --help` for a list
of all variables and options.

The following tools must be available:

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


# Getting started
PiGx consists of the following scripts:

| File          | Purpose       |
| ------------- |:-------------:|
| pigx-bsseq    |: (main script) - establishes config file, links to input, reference paths and launches. |
| [TableSheet].csv  |: (primary input file)  spreadsheet supplying basic parameters of the process: (e.g. filenames, paths, etc.) |
| BSseq_pipeline.py |  Defines the rules of the pipeline for data processing.     |
| config.json   | Generated automatically by pigx-bsseq. Defines various parameters; e.g. input/output folder paths, sample names, etc. |
| func_defs.py  | Subscript that defines various functions called in the main snakemake script.                          |

To run PiGx on your experimental data, first enter the necessary parameters in the spreadsheet file (see following section), and then from the terminal type

```
$ pigx-bsseq [options]
```

Upon doing so, the folder specified by "PATHOUT" in the input table will be created in addition to several subdirectories.
Most of these directories correspond to various stages of completion in the data-processing pipeline; however one additional folder will be called path_links/.
In this folder, symbolic links will be established that point directly to the reference genome and input files used for the calculation --provided they don't already exist. If they *do* already exist, then they will NOT be updated, and warnings to this effect will be suppressed. For this reason, you must BE CAREFUL WHEN RENAMING FILES IN THE INPUT FOLDER.
If, for example, you rename an input file to a name that has already been used (e.g. if you have "Sample1.fq.gz", "Sample2.fq.gz", and then delete Sample1 and decrement the names of subsequent files) or if you modify an extension (e.g. from ".fastq" to the standard ".fq") then the link can persist despite pointing to the wrong data.
If you absolutely must rename input files after having already run PiGx once, then go into ../PATHOUT/path_links/input/ and delete the corresponding links so that they can be generated from scratch automatically the next time PiGx is run. You have been warned.


# Input parameters

The input parameters specifying the desired behaviour of PiGx should
be entered into the tablesheet file.  When PiGx is run, the data from
this file will be used to automatically generate a configuration file
with the following values:
 
| Variable name | description |
| ------------- |:-----------:|
| NICE          | integer: from -20 to 19; higher values make the program execution less demanding on computational resources |
| bismark_args  | string: optional arguments supplied to bismark during alignment. See the [Bismark User Guide](https://rawgit.com/FelixKrueger/Bismark/master/Docs/Bismark_User_Guide.html#appendix-ii-bismark)  |
| SAMPLES       | struct: list of all the samples to be considered. |
| files         | string(s): part of SAMPLE: lists files (without extension) to read. when 2 are specified, paired-end is assumed, otherwise, single end. |
| PATHIN        | string: location of the experimental data files (.fastq[.gz\|.bz2])   |
| PATHOUT       | string: ultimate location of the output data and report files   |
| GENOMEPATH    | string: location of the reference genome data for alignment   |

 
All input files (paired or single end) must be present in the foler
indicated by `PATHIN`, And must have their files listed among
`SAMPLES`. Output from the snakemake script will then be sent to the
folder indicated by `PATHOUT`, with subdirectories corresponding to
the various stages of the process.

The directory indicated by `GENOMEPATH` must contain the reference
genome being mapped to.

