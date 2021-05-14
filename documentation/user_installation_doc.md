# Installation

This step by step installation and how-to-use guide should only allow test users to run the pipeline in this very preliminary state. The usage will change as soon the pipeline is more _PiGx-ifyied_.

## Download/Install the pigx-sarscov2-ww

Clone repository and enter the reproducible guix environment:

```
https://github.com/BIMSBbioinfo/pigx_sarscov2_ww.git
cd pigx_sarscov2_ww.git
USE_GUIX_INFERIOR=t guix environment -m manifest.scm
```

## Structure overview

This guide will lead to the following structure of directories and files.

```
tests/
|
|__sample_sheet.csv
|
|__settings.yaml
|
|__sample_data/__reads/
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
``` 


## Settings and Sample file

A draft for these two files are given in the `tests/` directory. The pre-filled rows should work but feel free to modify either locations in the settings file or given samples in the sample sheet. These files can not yet be provided as arguments in the final run command of the pipeline. Instead the `snakefile.py` relies on these two files to be in the `tests/` directory. 

## Dummy sample fastq file

For testing purposes we created small dummy reads `Test_R1.fastq` and `Test_R2.fastq` which can be found in `tests/sample_data/reads/`.

## Prepare databases

Before the pipeline can be run, 3 databases must be downloaded and their location will need to be provided in the settings file. Depending on the size of the databses this can take some time.
Be sure that the pigx-sarscov2-ww pipeline is downloaded and the tools are installed or used via the provided and suggested guix environment, bevore preparing the databases. One database (signature mutations, `sigmut_db`) is already provided via the repository. The folder structure is suggested like the following and pre-filled accordingly in the settings file.

### Kraken2 database

There are several libraries of genomes that can be used to classify the (unaligned) reads. It is up to you which one to use, but be sure that they fulfill the necessity stated by Kraken2 [Kraken2 manual](https://github.com/DerrickWood/kraken2/wiki/Manual#kraken-2-databases). We recommend to use the Plus-PFP library provided [here](https://benlangmead.github.io/aws-indexes/k2). 
It is also possible to have multiple Kraken2 databases, just be sure to provide the wanted one to the settings file.

First download and unpack the database in the `tests/databases/kraken_db/`:

```
mkdir tests/databases/kraken_db/
cd tests/databases/kraken_db/

# NOTE: This command will download a very large file, after unpacking this will
# require about 100GB of disc space. If this is not feasible use another
# database instead. For this please see link above or commented lines below.
 wget 'https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_20210127.tar.gz' --output-document 'k2_pluspfp_20210127.tar.gz'
 tar -xvf k2_pluspfp_20210127.tar.gz

# Use the following two lines to use smaller 8GB version instead.
# wget 'https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_8gb_20210127.tar.gz' --output-document 'k2_pluspfp_8gb_20210127.tar.gz'
# tar -xvf k2_pluspfp_8gb_20210127.tar.gz
```

Next move one folder up to `databases/` and build the Kraken database. This might take a while (depends on the size of the downloaded database):

```
cd ..
DBNAME=kraken_db/
kraken2-build --use-ftp --download-taxonomy --db $DBNAME # if this fails, you might want to try it without the --use-ftp flag
kraken2-build --build --db $DBNAME
```

Kraken might tell you that it can't find a library subdirectory. If that's the case it should be fine though. 

### Krona database

Krona Tools needs a two files, which have to be installed in the `tests/databases/krona_db/` folder. Also this might take a while:

```
DBNAME=krona_db/
mkdir $DBNAME
KRONA=$(dirname $(which ktImportTaxonomy))/../share/krona-tools/ # this is just a workaround until tool paths are declared
$KRONA/updateTaxonomy.sh $DBNAME # the scripts are stored a priori in that folder
$KRONA/updateAccessions.sh $DBNAME
```

### VEP database

Just download sars_cov_2 database for VEP and unpack it in the `tests/databases/vep_db/` directory.

```
DBNAME=vep_db/
mkdir $DBNAME
cd $DBNAME
wget ftp://ftp.ensemblgenomes.org/pub/viruses/variation/indexed_vep_cache/sars_cov_2_vep_101_ASM985889v3.tar.gz --output-document 'sars_cov_2_vep_101_ASM985889v3.tar.gz'
tar -xvf sars_cov_2_vep_101_ASM985889v3.tar.gz
```

### sigmut database

Nothing to be done here. Necessary files are provided in `tests/databases/sigmut_db/`.

## run pipeline

As mentioned above this how-to-use guide will change in the near future. This step in particular. Right now it is done by manually running snakemake while providing the `snakefile.py` and at least the number of jobs. 

```
NUM_JOBS = 4
snakemake --snakefile=snakefile.py --jobs=$NUM_JOBS
```
