<a name="logo"/>
<div align="center">
<img src="images/Logo_PiGx.png" alt="PiGx Logo"  width="30%" height="30%" ></img>
</a>
</div>

**Copyright 2017-2022: Vedran Franke, Alexander Blume, Ricardo Wurmus.**
**This work is distributed under the terms of the GNU General Public License, version 3 or later.  It is free to use for all purposes.**

-----------

## Summary

PiGx ChIPseq (pipelines in genomics for Chromatin Immunoprecipitation
Sequencing) is an analysis pipeline for preprocessing, peak calling and
reporting for ChIP or ATAC sequencing experiments. It is easy to use and produces high
quality reports. The inputs are reads files from the sequencing experiment, and
a configuration file that describes the experiment. In addition to quality
control of the experiment, the pipeline enables multiple peak calling
analysis and allows the generation of a UCSC track hub in an easily
configurable manner. 
Starting from *versions >= 0.0.50* the pipeline allows to perform 
highly customizable differential analysis using DEDeq2 to detect genome regions with 
differences in occupancy between defined groups.



## What does it do

- Trim reads using trim-galore 
- Quality control reads using fastQC and multiQC
- Map reads to genome using Bowtie2
- Call peaks for multiple combinations of samples using MACS2
- Control reproducibility of experiments using IDR
- Generate a UCSC track hub to view in Genome Browser  
- *NEW* Highly customizable Differential Analysis


## What does it output

- QC reports
- bam files
- bigwig files
- narrowPeak files
- UCSC track hub folder
- *NEW* Differential Analysis Report

# Installation of the PiGx-ChIPseq Pipeline

## Install via guix

You can install this pipeline with all its dependencies using GNU Guix:

    guix package -i pigx-chipseq
    
After exporting the guix profile to the PATH you should be able to [get started](#getting-started).

## Install from source

You can also install it from source manually.  You can find the
[latest
release](https://github.com/BIMSBbioinfo/pigx_chipseq/releases/latest)
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
Guix](https://gnu.org/s/guix). If you do not use one of these package
managers, you will need to ensure that the following software is
installed:

<details>
<summary>Software dependencies</summary>

- R
    - argparser
    - biocparallel
    - biostrings
    - chipseq
    - data.table
    - dyplr
    - genomation
    - genomicalignments
    - genomicranges
    - rsamtools
    - rtracklayer
    - rsubread
    - s4vectors
    - stringr
    - jsonlite
    - heatmaply
    - htmlwidgets
    - ggplot2
    - ggrepel
    - plotly
    - rmarkdown
    - ggrepel
    - deseq2
    - dt
    - pheatmap
    - corrplot
    - reshape2
    - scales
    - crosstalk
    - gprofiler2
    - summarizedexperiment
    - hexbin
- python
    - snakemake
    - pyyaml
    - wrapper
    - pytest
    - xlrd
    - magic
- pandoc
- fastqc
- multiqc
- trim-galore
- bowtie
- macs2
- idr
- samtools
- bedtools
- bedToBigBed
- bamToBed


</details>


## via Guix

Assuming you have Guix installed, the following command spawns a
sub-shell in which all dependencies are available:

```sh
guix environment -l guix.scm
```


# Getting started

To run PiGx on your experimental data, first enter the necessary parameters in the spreadsheet file (see following section), and then from the terminal type

```sh
$ pigx-chipseq [options] sample_sheet.csv
```
To see all available options type the `--help` option

```sh
$ pigx-chipseq --help

usage: pigx-chipseq [-h] [-v] -s SETTINGS [-c CONFIGFILE] [--target TARGET]
                   [-n] [--graph GRAPH] [--force] [--reason] [--unlock]
                   samplesheet

PiGx ChIPseq Pipeline.

PiGx ChIPseq is a data processing pipeline for ChIPseq read data.

positional arguments:
  samplesheet                             The sample sheet containing sample data in yaml format.

optional arguments:
  -h, --help                              show this help message and exit
  -v, --version                           show program version number and exit
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

## Input Parameters

The pipeline requires two files as input to specify the [samples](#sample-sheet)  and the design of the [analysis](#settings-file).  

### Sample Sheet

The samples used for any subsequent analysis are defined in the _sample sheet_ section. 

| SampleName | Read | Read2 |
|------|-------|--------|

- SampleName is the name for the sample
- _Read/Read2_ are the fastq file names of paired end reads
  - the location of these files is specified in `settings.yaml`
  - for single-end data, leave the Read2 column in place, but have it empty

#### Technical Replicates

The sample sheet offers support for technical replicates, by repeating the sample name (first column) for different input files (second,third column).
The quality check will be performed for any input file and replicates will be merged during the mapping. 


| SampleName | Read               | Read2              | *Group (optional)* |
| ------     | -------            | --------           | --------           |
| ChIPpe     | ChIPpe_R1.fq.gz    | ChIPpe_R2.fq.gz    | Group1             |
| ChIPpe     | ChIPpe_t2_R1.fq.gz | ChIPpe_t2_R2.fq.gz | Group2             |

### Settings File

The settings file is a file in yaml format specifying general settings and the details of the analysis. It has the following **required** sections: 

#### Locations

Defines paths to be used in the pipeline, some of the items are required and some optional (can stay blank): 

| item    | required | description |
|---------|----------|-------------|
| _input-dir_ | yes | directory of the input files (`fastq` files) |
| _output-dir_    | yes | output directory for the pipeline |
| _genome-file_    | yes | path to the reference genome in `fasta` forma |
| _index-dir_    | no | directory containing pre-built mapping indices for the  given reference genome (created with `bowtie2-build`) |
| _gff-file_    | no | location of a `GTF` file with genome annotations for the  given reference genome |

#### General

These are settings which apply to all analysis (unless adjusted in single analysis):

| item    | required | description |
|---------|----------|-------------|
| _assembly_ | yes | version of reference genome (e.g. hg19,mm9, ...) |
| *spikein_name* | no | version of spike-in genome (e.g. hg19,mm9, ...) |
| *organism* | no | needed for GO term analysis (e.g. hsapiens, mmusculus, dmelanogaster, celegans) |
| _params_    | no | list of default parameters for tools and scripts (for tools check respective manual for available parameters) |

#### Execution

The `execution` section in the settings file allows the user to specify whether the pipeline is to be submitted to a cluster, or run locally, and the degree of parallelism. For a full list of possible parameters, see `etc/settings.yaml`.


A minimal settings file could look like this, but please consider that no analysis will be performed without adding [analysis information](#analysis-sections) :

```yaml
locations:
  input-dir: in/reads/
  output-dir: out/
  genome-file: genome/my_genome.fa
  index-dir:
  gff-file: genome/mm_chr19.gtf

general:
  assembly: mm9
  spikein_name: hg38
  # The "organism" field is needed for GO term analysis. Leave it empty
  # if not interested in GO analysis.  Otherwise provide a string with
  # the initial of genus and the species name (e.g. hsapiens, mmusculus,
  # dmelanogaster, celegans)
  organism: 'mmusculus'
  params:
    export_bigwig:
        extend: 200
        scale_bw: 'yes'
    bowtie2:
        N: 0
        k: 1
    idr:
        idr-threshold: 0.1
    macs2:
        #        g: mm
        keep-dup: auto
        q: 0.05
        nomodel: ''
    extract_signal:
        expand_peak: 200
        bin_num: 20
    feature_counting:
        minFragLength: 50
        maxFragLength: 600
        countChimericFragments: ''
        requireBothEndsMapped: ''

execution:
  submit-to-cluster: no
  rules:
    __default__:
      queue: all.q
      memory: 8G
    bowtie2:
      queue: all.q
      memory: 16G

```

### Analysis Sections

The analysis part of the setting file describes the experiment. It has following sections: 

| section | required | description |
|---------|----------|-------------|
| *peak_calling*  | yes | defines which samples will be used to detect regions of enriched binding ( multiple combinations and variations are possible, [see here for details](#peak-calling) ) |
| _idr_ | no | specifies pairs of *peak calling* analysis that are compared to determine the reproducibilty of the general experiment ([see here for details](#optionalidr)) |
| _hub_ | no | describes the general layout of a UCSC hub that can be created from the processed data and allows the visual inspection of results at a UCSC genome browser ([see here for details](#optional-hub)) |
| *feature_combination* | no | defines for a list of *peak calling* and/or *idr* analysis the combination of regions shared among this list ([see here for details](#optional-feature-combination)) |
| *differential_analysis* | no | specifies for a list of groups defined in the sample sheet which samples are used as Case or Control in the Differential Analysis  ([see here for details](#optional-differentialr-analysis)) |


The creation of these sections is straight forward considering the following snippets as template. Comments and examples within the snippets provide guidance of what is possible and what to take care of.

#### Peak Calling

The previously defined samples are used for subsequent peak calling analysis to detect regions of enriched binding. In this section any number of comparisons can be defined, while multiple combinations and variations are possible. In terms of peak calling the **ChIP** (also called treatment) is the sample in which we want to detect enriched regions compared to the **Cont(rol)** (or background) sample. Each analysis can be run with a unique set of parameters and default parameters for all analysis can be defined in the [settings file](#settings-file) , check available parameters and description [here](https://github.com/taoliu/MACS).
For more information have a look at the publication for the software we are using "Zhang et al. Model-based Analysis of ChIP-Seq (MACS). Genome Biol (2008) vol. 9 (9) pp. R137".  

```yaml
# define peak calling analysis
peak_calling:
    # analysis can have any name, but the names have to be unique 
    Peaks1: 
        # sample(s) to be used as treatment sample 
        ChIP: ChIP1
        # sample(s) to be used as control sample
        Cont: Cont1
        params:
            macs2:
                # each analysis can be adjusted independently
                # add/modify available parameters of the analysis
                nomodel: ''
                extsize: 300
    Peaks2:
        ChIP: ChIP2
        Cont: Cont2
        params:
            macs2:
                # each analysis can be adjusted independently
                nomodel: ''
                extsize: 147
    Peaks4:
        ChIP:
            # multiple samples can be used as treatment
            - ChIP1
            - ChIP2
        Cont:
            # multiple samples can be used as control
            - Cont1
            - Cont2
        params:
            macs2:
                nomodel: ''

    Peaks5:
        # the number of samples per group can differ 
        ChIP: ChIP2
        Cont:
            - Cont1
            - Cont2
        params:
            macs2:
                nomodel: ''

    Peaks6:
        # analysis can be performed without control
        ChIP: ChIP1
        Cont:
        params:
            macs2:
                nomodel: ''
```

#### (_optional_) IDR

Assuming that the some samples are (biological/technical) replicates, in order to measure the consistency between them use the irreproducible discovery rate (IDR)  "Li, Q., Brown, J. B., Huang, H., & Bickel, P. J. (2011). Measuring reproducibility of high-throughput experiments. The annals of applied statistics, 5(3), 1752-1779.", which is in general a good (but very stringent) quality control.

```yaml
idr:
    # idr analysis can have any name, but the names have to be unique 
    ChIP_IDR:
        # define the pair of samples, add more combinations for more replicates
        # ChIP1 and ChIP2 are required labels, only adjust values of those
        ChIP1: Peaks1
        ChIP2: Peaks2
```

#### (_optional_) Hub

In the _hub_ section the general layout of a [UCSC Track Hubs](https://genome.ucsc.edu/goldenpath/help/hgTrackHubHelp.html#Intro) is described with some minimal items. The track hub is generated from the processed data and allows the visual inspection of results at a UCSC genome browser (for supported genomes). 

The required items to define the hub are the following:

| item    | example | description |
|---------|----------|-------------|
| name    | PiGx_Hub | name of the hub directory |
| shortLabel    | PiGx_Short | short name of hub is displayed as name above track groups |
| longLabel    | PiGx_Hub_Long | descriptive longer label for hub is displayed as hub description |
| email    | my.mail[at]domain.com | whom to contact for questions about the hub or data |
| descriptionUrl    | pigx_hub.html | URL to HTML page with a description of the hub's contents |
| super_tracks    | see below | specification of hub layout (track groups, tracks) |

This is a small example how this could look like:

```yaml
hub:
    name: PiGx_Hub
    shortLabel: PiGx_Short
    longLabel: PiGx_Hub_Long
    email: my.mail@domain.com
    descriptionUrl: pigx_hub.html
    super_tracks:
        # track groups can have any name, but the names have to be unique 
        Tracks1:
            # tracks can have any name, but the names have to be unique 
            track11:
                # to add peaks as a track, define "type: macs" 
                name: Peaks1
                type: macs
            track12:
                # to add coverage signal as a track, define "type: bigwig"
                name: ChIP1
                type: bigWig
            # descriptive longer label for track group is
            # displayed as description in track settings
            long_label: Tracks1_long
```

#### (_optional_) Feature Combination

To find the combination of enriched binding regions, which is shared among a set of *peak calling* and/or *idr* analysis results, define a feature in the *feature_combination* section. Only items defined in the *peak_calling* and *idr* sections can be used here.  

```yaml
feature_combination:
    # features can have any name, but the names have to be unique
    Feature1:
        # define feature based on only one result
        - ChIP_IDR
    Feature2:
        # define feature based on more than one result
        - Peaks6
        - Peaks5
    Feature3:
        # define feature based on different analysis types
        - ChIP_IDR
        - Peaks5
```

#### (_optional_) Differential Analysis

To find genomic regions that are differentially bound by a certain chipped protein or that show differential open chomatin patterns, define a differential analysis to compare *Case* and *Control* samples, with optional covariates, using the DESeq2 framework ("Love, M.I., Huber, W. and Anders, S., 2014. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome biology, 15(12), p.550."). Optionally predefined *peak calling* and/or *idr* analysis resultscan be used to specify pre-called features that are combined to regions of interest, but by default a *peak calling* analysis containing all samples of the specifed groups is performed to create the consensus peaks defining the regions of interest. Each analysis can be run with a unique set of parameters for the peak calling and read summarization, however default parameters for all analysis can be defined in the general section of the [settings file](#settings-file). Reads over features are counted using Rsubreads featureCount function, as this allows to finetune the summarization by multiple options explained in the [package's manual](https://www.bioconductor.org/packages/release/bioc/vignettes/Rsubread/inst/doc/SubreadUsersGuide.pdf), please see their publication for more details: "Liao, Y., Smyth, G.K. and Shi, W., 2019. The R package Rsubread is easier, faster, cheaper and better for alignment and quantification of RNA sequencing reads. Nucleic acids research, 47(8), pp.e47-e47.". 

```yaml
differential_analysis:
    Analysis1:
        Case:
            - Group1 
            - Group2
        Control:
            - Group3
        # comma separated list of additional co-variates to control for
        # in differential expression analysis (e.g. batch, age,
        # temperature, sequencing_technology etc.). Must correspond to a
        # column field in the sample_sheet.csv file.)    
        Covariates: ''
        Peakset: # choose peaks to define consensus peak set
            - ChIP_IDR
            - Peaks6
            - Peaks5
        params:
            feature_counting:
                minFragLength: 50
                maxFragLength: 150
    Analysis2:
        Case:
            - Group1 
            - Group2
        Control:
            - Group3
        Peakset: # if no peaks given call joint peaks for given samples
    Analysis3:
        Case: 
            - Group1 
            - Group2
        Control: Group3
        params:
            macs2:
                nomodel: ''
                extsize: 147
```

## Output Folder Structure

```
|-- Analysis
|-- Annotation
|-- BigWig
|-- Bowtie2_Index
|-- FastQC
|-- Log
|-- Mapped
|-- Peaks
|-- Reports
|-- Trimmed
|-- UCSC_HUB
```

#### Analysis

Contains RDS files with intermediary analysis steps. RDS are binary files which efficiently 
store R objects.

#### Annotation

Formatted GTF annotation.

#### BigWig

Symbolic links to the bigWig signal files.

#### Bowtie2_Index

Processed genme file along with the Bowtie2_Index

#### FastQC

FastQC sequencing quality report

#### Log 

Detailed output from execution of each step of the pipeline.

#### Mapped

Mapped reads in .bam format, and corresponding bigWig files.

#### Peaks

Peaks called with MACS2. Depending on the parameters, contains either narrowPeak or 
broadPeak format. **sample_qsort.bed** contains uniformly processed peaks, sorted by
their corresponding p value.

#### Reports

Contains MultiQC and ChIP quality reports in html format.
Optionally contains folders with Results of Differential Analysis  wth 
DESeq2 and featureCounts tables in CSV and Reports in html format.

#### Trimmed

Trimgalore adaptor and quality trimmed files.

#### UCSC_Hub

Contains a completely formatted UCSC hub, with track descriptions, peaks and bigWig tracks.

# Documentation

If are interested in the general workflow of the pipeline and want to learn more about PiGx make sure you check:
- the docs at https://bioinformatics.mdc-berlin.de/pigx_docs/pigx-chip-seq.html#pigx-chip-seq
- our other pipelines at https://bioinformatics.mdc-berlin.de/pigx/

# Citation

If you use the pipeline in published research, please cite our publication:

Wurmus, R., Uyar, B., Osberg, B., Franke, V., Gosdschan, A., Wreczycka, K., Ronen, J., Akalin, A. (2018). PiGx: Reproducible genomics analysis pipelines with GNU Guix. GigaScience. https://doi.org/10.1093/gigascience/giy123 

# Questions
If you have further questions please e-mail:
pigx@googlegroups.com or use the web form to ask questions
https://groups.google.com/forum/#!forum/pigx/
