# pigx_chipseq - enables easy and reproducible analysis of ChIP-seq data

pigx_chipseq is a snakemake pipeline wrapped in GUIX, which enables
easy and reproducible analsys of ChIP-seq data.

## Installation

## Usage

## Input Parameters

### Sample Sheet

The sample sheet is a file in yaml format describing the experiment. It has multiple sections: 

| section | required | description |
|---------|----------|-------------|
| _samples_ | yes | describes the mapping of samples, specifying read file names and library type (_single_/_paired_) ([see here for details](#samples)) |
| *peak_calling*  | yes | defines which samples will be used to detect regions of enriched binding ( multiple combinations and variations are possible, [see here for details](#peak_calling) ) |
| _idr_ | no | specifies pairs of *peak calling* analysis that are compared to determine the reproducibilty of the general experiment ([see here for details](#idr)) |
| _hub_ | no | describes the general layout of a UCSC hub that can be created from the processed data and allows the visual inspection of results at a UCSC genome browser ([see here for details](#hub)) |
| *feature_combination* | no | defines for a list of *peak calling* and/or *idr* analysis the combination of regions shared among this list ([see here for details](#feature_combination)) |


#### Sample Sheet configuration

##### Samples

The samples used for any subsequent analysis are defined in the _samples_ section. 

```
# define mapping
samples:
    # samples can have any name, but the names have to be unique
    ChIP1: 
        fastq:
            # file names of raw data in fastq format, either gzipped or not 
            - ChIP.fq.gz
        # map reads in single-end mode
        library: single
    Cont1:
        fastq:
            - Cont.fq.gz
        library: single
    ChIPpe:
        fastq:
            - ChIPpe_R1.fq.gz
            - ChIPpe_R2.fq.gz
        # map reads in paired-end mode
        library: paired
```

##### Peak Calling

The previously defined samples are used for subsequent peak calling analysis to detect regions of enriched binding. In this section any number of comparisons can be defined, while multiple combinations and variations are possible. In terms of peak calling the **ChIP** (also called treatment) is the sample in which we want to detect enriched regions compared to the **Cont(rol)** (or background) sample. Each analysis can be run with a unique set of parameters and default parameters for all analysis can be defined in the [settings file](#settings-file) , check available parameters and description [here](https://github.com/taoliu/MACS).
For more information have a look at the publication for the software we are using "Zhang et al. Model-based Analysis of ChIP-Seq (MACS). Genome Biol (2008) vol. 9 (9) pp. R137".  

```
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

##### (_optional_) IDR

Assuming that the some samples are (biological/technical) replicates, in order to measure the consistency between them use the irreproducible discovery rate (IDR)  "Li, Q., Brown, J. B., Huang, H., & Bickel, P. J. (2011). Measuring reproducibility of high-throughput experiments. The annals of applied statistics, 5(3), 1752-1779.", which is in general a good (but very stringent) quality control.

```
idr:
    # idr analysis can have any name, but the names have to be unique 
    ChIP_IDR:
        # define the pair of samples, add more combinations for more replicates
        ChIP1: Peaks1
        ChIP2: Peaks2
```

##### (_optional_) Hub

In the _hub_ section the general layout of a [UCSC Track Hubs](https://genome.ucsc.edu/goldenpath/help/hgTrackHubHelp.html#Intro) is described with some minimal items. The track hub is generated from the processed data and allows the visual inspection of results at a UCSC genome browser (for supported genomes). 

The required items to define the hub are the following:

| item    | example | description |
|---------|----------|-------------|
| name    | Pix_Hub | name of the hub directory |
| shortLabel    | Pix_Short | short name of hub is displayed as name above track groups |
| longLabel    | Pix_Hub_Long | descriptive longer label for hub is displayed as hub description |
| email    | my.mail[at]domain.com | whom to contact for questions about the hub or data |
| descriptionUrl    | pix_hub.html | URL to HTML page with a description of the hub's contents |
| super_tracks    | see below | specification of hub layout (track groups, tracks) |

This is a small example how this could look like:

```
hub:
    name: Pix_Hub
    shortLabel: Pix_Short
    longLabel: Pix_Hub_Long
    email: vedran.franke@mdc-berlin.de
    descriptionUrl: pix_hub.html
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

##### (_optional_) Feature Combination

To find the combination of enriched binding regions, which is shared among a set of *peak calling* and/or *idr* analysis results, define a feature in the *feature_combination* section. Only items defined in the *peak_calling* and *idr* sections can be used here.  

```
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

### Settings File

The settings file is a file in yaml format specifying general settings. It has the following sections: 

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
| _params_    | no | list of default parameters for tools and scripts (for tools check respective manual for available parameters) |

#### Execution

The `execution` section in the settings file allows the user to specify whether the pipeline is to be submitted to a cluster, or run locally, and the degree of parallelism. For a full list of possible parameters, see `etc/settings.yaml`.


The settings.file will look like this:

```
locations:
  input-dir: in/reads/
  output-dir: out/
  genome-file: genome/my_genome.fa
  index-dir:
  gff-file: genome/mm_chr19.gtf

general:
  assembly: hg19
  params:
    extend: 200
    scale_bw: 'yes'
    bowtie2:
        k: 1
    idr:
        idr-threshold: 0.1
    macs2:
        g: hs
        keep-dup: auto
        q: 0.05
    extract_signal:
        expand_peak: 200
        bin_num: 20

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





## Reports


