# Introduction 
PiGx SARS-CoV-2 is a pipeline for analysing data from sequenced wastewater samples and identifying given variants-of-concern of SARS-CoV-2. Currently wastewater samples are used, which are enriched for SARS-CoV-2. The pipeline can be used for continuous sampling. The output of the PiGx SARS-CoV-2 pipeline is summarized in a report which provides an intuitive visual overview about the development of variant abundance over time and location. Additionally there will be more detailed reports per sample, which cover the quality control of the samples, the detected variants and a taxonomic classification of all reads which are not aligned to SARS-CoV-2.

## Workflow

First the raw reads are trimmed by using [Prinseq](http://prinseq.sourceforge.net/) to improve alignment rates and mutation calling. Next, the trimmed reads are aligned to the reference genome of SARS-CoV-2 using [BWA](https://github.com/lh3/bwa), and the results are *SAM*/*BAM* files of **aligned** and **unaligned reads**. Following the alignment a quality check on raw and processed reads is performed by using [MultiQC](https://multiqc.info/). Calling the variants and inferring SNVs (single nucleotide polymorphisms) on the **aligned reads** is done with [LoFreg](https://csb5.github.io/lofreq/). Mutations are annotated with [VEP](https://covid-19.ensembl.org/index.html). Estimation of Variants of Concern (VOC) frequencies are done by deconvolution. 
To investigate the abundance of other existing species in the wastewater samples the **unaligned reads** will be taxonomicly classified with [Kraken2](https://github.com/DerrickWood/kraken2). The Kraken2 requires a database of the genomes against the reads are getting aligned, therefore keep in mind that you can only find those species which are included in the chosen database. For documentation how to set this up, see: [Prepare databases](#prepare-databases). For a better and interactive visualization of all species present in the wastewater [Krona](https://github.com/marbl/Krona/wiki) is used. Also here a small step of setting up a database is needed before running the pipeline, see: [Prepare databases](#prepare-databases). 

## Output
* Merged Report including:
   * Overview of development of variant and mutation abundance over time and locations
   * Quality Control report of raw and processed (trimmed) reads
   * Variant report
   * Taxonomic classification
* SAM/ BAM files of the aligned and unaligned reads against SARS-CoV-2


# Installation

## Install via guix

This pipeline is going to be packaged in the [GNU Guix](https://gnu.org/s/guix) package manager soon. Right now in this preliminary state it has to be installed from source manually.

## Install from source

First, please clone this repository and change directory accordingly:

```sh
git clone https://github.com/BIMSBbioinfo/pigx_sarscov2_ww.git
cd pigx_sarscov2_ww
```

To fetch code that is common to all PiGx pipelines run this:

```sh
git submodule update --init
```

Before setting everything up, though, make sure all dependencies are met by either installing the following software manually, or by entering the provided reproducible Guix environment. If you are using Guix we definitely recommend the latter. This command spawns a sub-shell in which all dependencies are available:

```sh
USE_GUIX_INFERIOR=t guix environment --pure -m manifest.scm --preserve=GUIX_LOCPATH
```

To use your current Guix channels instead of the fixed set of
channels, just omit the `USE_GUIX_INFERIOR` shell variable:

```sh
guix environment --pure -m manifest.scm --preserve=GUIX_LOCPATH
```

Note that `--pure` unsets all environment variables that are not
explicitly preserved.  To access other executables that are not part
of the environment please address them by their absolute file name.


<details>
<summary>Software dependencies</summary>

- R
    - minimal
    - base64url
    - dplyr
    - dt
    - ggplot2
    - magrittr
    - plotly
    - qpcr
    - rmarkdown
    - stringr
    - tidyr
    - reshape2
- python
    - wrapper
    - pyyaml
- bash-minimal
- bwa
- ensembl-vep
- fastqc
- multiqc
- kraken2
- krona-tools
- lofreq
- prinseq
- samtools
- snakemake
- autoconf
- automake
- coreutils
- gawk
- grep
- make
- sed


</details>
</br>

Inside the environment you can then perform the usual build steps:

```sh
./bootstrap.sh # to generate the "configure" script
./configure
make
make check
```

# Getting started

At this point you are able to run the pipeline from within the current directory `pigx_sarscov2_ww`. Use `--help` to see the available options.  To ensure that PiGx uses files from the source directory set the shell variable `PIGX_UNINSTALLED` to any value.

```sh
PIGX_UNINSTALLED=t ./pigx-sars-cov2-ww --help
```
<details>
<summary>toggle output</summary>

```sh
usage: pigx-sars-cov2-ww [-h] [-v] (--init [{settings,sample-sheet,both}] | -s SETTINGS) [-c CONFIGFILE]
                         [--target TARGET] [-n] [--graph GRAPH] [--force] [--reason] [--unlock] [--verbose]
                         [--printshellcmds]
                         [sample_sheet]

PiGx SARS-CoV-2 wastewater sequencing Pipeline.

This is a pipeline for analyzing data from sequenced wastewater
samples and identifying given variants-of-concern of SARS-CoV-2.  The
pipeline can be used for continuous sampling.  The output report
provides an intuitive visual overview about the development of variant
abundance over time and location.

positional arguments:
  sample_sheet                            The sample sheet containing sample data in CSV format.

optional arguments:
  -h, --help                              show this help message and exit
  -v, --version                           show program's version number and exit
  --init [{settings,sample-sheet,both}]   Generate a template SETTINGS file, a SAMPLE-SHEET.  Leave empty for both.
  -s SETTINGS, --settings SETTINGS        A YAML file for settings that deviate from the defaults.
  -c CONFIGFILE, --configfile CONFIGFILE  The config file used for calling the underlying snakemake process.  By
                                          default the file 'config.json' is dynamically created from the sample
                                          sheet and the settings file.
  --target TARGET                         Stop when the named target is completed instead of running the whole
                                          pipeline.  The default target is "final-report".  Pass "--target=help"
                                          to describe all available targets.
  -n, --dry-run                           Only show what work would be performed.  Do not actually run the
                                          pipeline.
  --graph GRAPH                           Output a graph in PDF format showing the relations between rules of
                                          this pipeline.  You must specify a graph file name such as
                                          "graph.pdf".
  --force                                 Force the execution of rules, even though the outputs are considered
                                          fresh.
  --reason                                Print the reason why a rule is executed.
  --unlock                                Recover after a snakemake crash.
  --verbose                               Print supplementary info on job execution.
  --printshellcmds                        Print commands being executed by snakemake.

This pipeline was developed by the Akalin group at MDC in Berlin in 2017-2021.
```
</details>
</br>

Though, to actually use it on your experimental data still more setup is required. Please follow [the steps to prepare the required databases](#prepare-databases) first and [prepare the input](#preparing-the-input). Then afterwards, you can run the pipeline from the source directory.

```sh
PIGX_UNINSTALLED=t ./pigx-sars-cov2-ww [options] sample_sheet.csv
```

For example, with this command the pipeline is used for the included test data: 

```sh
PIGX_UNINSTALLED=t ./pigx-sars-cov2-ww -s tests/settings.yaml tests/sample_sheet.csv
```

## Preparing the input

In order to run the pipeline, the user must supply

- a sample sheet
- a settings file

both files are described below.

In order to generate template settings and sample sheet files, type

```sh
PIGX_UNINSTALLED=t ./pigx-sars-cov2-ww --init
```

in the shell, and a boilerplate `sample_sheet.csv` and `settings.yaml` will be written to your current directory. An example for both files is provided in the `tests/` directory.

### Sample sheet

The sample sheet is a tabular file (`csv` format) describing the experiment. The table has the following columns:

| SampleName | Read               | Read2              | date               | location_name      | coordinates_lat    | coordinates_long   |
| ------     | -------            | --------           | --------           | --------           |  --------          |  --------          |
| Test0      | Test0_R1.fastq     | Test0_R2.fastq     | 2021-01-01T08:00:00| Berlin             | 52.3646054650197   | 13.5098643274129   |
| Test2      | Test2_R1.fastq     | Test2_R2.fastq     | 2021-01-03T08:00:00| Munich             | 48.2084486780314   | 11.6282300407931   |

- _SampleName_ is the name for the sample
- _Read_ & _Read2_ are the fastq file names of paired end reads
  - the location of these files is specified in the settings file
  - single-end data is not yet supported
  - compressed (`.gz`) reads are not yet supported
- _date_ is a date/time in ISO format (`yyyy-mm-ddThh:mm:ss`)
- _location_name_ is the name of the location and should be unique per coordinates
- _coordinates_lat_ & _coordinates_long_ correspond the latitude and longitude of the location name

### Settings file

In the settings file, parameters are saved, in YAML format, to configure the execution of PiGx-sarscov2. It specifies:

#### Locations

- _output-dir_, the location of the outputs for the pipeline
- _reads-dir_, the location of the reads (directory where `fastq` files are)
- _reference-fasta_, the `fasta` file with the reference genome (must be prepared by the user)
- _amplicons-bed_, the amplicons `bed` file for coronavirus (must be prepared by the user)
- _kraken-db-dir_, the location of the kraken database (must be prepared by the user)
- _krona-db-dir_, the location of the krona database (must be prepared by the user)
- _sigmut-db-dir_, the location of the signature mutations database (provided at databases/sigmut_db/)
- _vep-db-dir_, the location of `sars_cov_2` database for VEP (must be prepared by the user)

## Prepare databases

Before the pipeline can work, three databases must be downloaded and their location will need to be provided in the settings file. Depending on the size of the databases this can take some time.
Be sure that the pigx-sarscov2-ww pipeline is downloaded and the tools are installed or used via the provided and suggested guix environment. One database (signature mutations, `sigmut_db`) is already provided via the repository. The directory structure is suggested like [this](#structure-overview) and pre-filled accordingly in the settings file. 

### Kraken2 database

There are several libraries of genomes that can be used to classify the (unaligned) reads. It is up to you which one to use, but be sure that they fulfill the necessities stated by Kraken2 [Kraken2 manual](https://github.com/DerrickWood/kraken2/wiki/Manual#kraken-2-databases). For an overall overview we recommend to use the Plus-PFP library provided [here](https://benlangmead.github.io/aws-indexes/k2). If the classification is not of concern or only the viruses are of interest, we recommend to use a smaller one. This will accelerate the speed.
It is also possible to have multiple Kraken2 databases installed, just be sure to provide the correct location to the settings file.

First download and unpack the database in the `databases/kraken_db/`:

```
DIR=databases/kraken_db/
mkdir -p $DIR

# NOTE: This command will download a very large file, after unpacking this will
# require about 100GB of disc space. If this is not feasible use another
# database instead. For this please see link above or commented lines below.
wget -qO- https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_20210127.tar.gz | tar -C $DIR -xzv

# Use the following two lines to use smaller 8GB version instead.
# wget -qO- https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_8gb_20210127.tar.gz | tar -C $DIR -xzv
```

Next go to the database/ directory and build the Kraken database. This might take a while (depends on the size of the downloaded database):

```
#current location: pigx_sarscov2_ww/
cd database/
DBNAME=kraken_db
kraken2-build --use-ftp --download-taxonomy --db $DBNAME # if this fails, you might want to try it without the --use-ftp flag
kraken2-build --build --db $DBNAME
```

Kraken might tell you that it can't find a library subdirectory. If that's the case it should be fine though.  
This folder should now contain at least these files: hash.k2d, opts.k2d and taxo.k2d.

### Krona database

Krona Tools needs two files, which have to be installed in the `databases/krona_db/`. Also this might take a while:

```
#current location: pigx_sarscov2_ww/databases/
DBNAME=krona_db/
mkdir -p $DBNAME
KRONA=$(dirname $(which ktImportTaxonomy))/../share/krona-tools/ # this is just a workaround until tool paths are declared
$KRONA/updateTaxonomy.sh $DBNAME # the scripts are stored a priori in that folder
$KRONA/updateAccessions.sh $DBNAME
```

### VEP database

Just download the `SARS_CoV_2` database for VEP (variant effect predictor) and unpack it in the `databases/vep_db/` directory.

```
#current location: pigx_sarscov2_ww/databases/
DBNAME=vep_db/
mkdir -p $DBNAME
wget -qO- ftp://ftp.ensemblgenomes.org/pub/viruses/variation/indexed_vep_cache/sars_cov_2_vep_101_ASM985889v3.tar.gz | tar -C $DBNAME -xzv
```

### Sigmut database

Necessary files are provided in `databases/sigmut_db/` for the current main Variants of Concern. Users can add files with new variants if/when necessary.


# Quick Start

To check wether the pipeline together with the databases was properly installed, run PiGx SARS-CoV-2 Wastewater Sequencing Pipeline on a minimal test dataset.
For this there are samples provided in `tests/sample_data/`. The folder structure should be provided like this, assuming all databases are set up like described [here](#prepare-databases):


```
pigx_sarscov2_ww
|
|__databases/
   │
   │__kraken_db/
   │
   │__krona_db/
   │
   │__vep_db/
   |
   |__sigmut_db/
|      
|__tests/
   |__sample_sheet.csv
   |
   |__settings.yaml
   |
   |__sample_data/
      |
      |__reads/
      |
      |__NC_045512.2.fasta
      |
      |__nCoV-2019_NCref.bed
|
|__[...]

``` 

Now the the test set can be run with the command: 

```
#current location: pigx_sarscov2_ww/
PIGX_UNINSTALLED=t ./pigx-sars-cov2-ww -s tests/settings.yaml tests/sample_sheet.csv
```

Inside `tests/` a new directory `output` is created, which includes specific directories containing output data for the respective step of the pipeline.    The `tests/output/reports/index.html` gives the overview over all merged reports for the test data. 

# Troubleshooting

If you have any questions please e-mail: pigx@googlegroups.com or use the web form to ask questions https://groups.google.com/forum/#!forum/pigx/. 

If you run into any bugs, please open an issue here: https://github.com/BIMSBbioinfo/pigx_rnaseq/issues. 

