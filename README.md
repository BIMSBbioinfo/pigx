<a name="logo"/>
<div align="center">
<img src="images/Logo_PiGx.png" alt="PiGx Logo"  width="30%" height="30%" ></img>
</a>
</div>

# PiGx SARS-CoV-2 Wastewater Sequencing Pipeline

**Copyright 2021-2022: Vic-Fabienne Schumann, Ricardo Wurmus, Miriam Faxel, Jan Dohmen, Rafael Cuadrat, Bora Uyar, Vedran Franke, Alexander Blume, and Altuna Akalin.**
**This work is distributed under the terms of the GNU General Public License, version 3 or later.  It is free to use for all purposes.**

-----------

PiGx SARS-CoV-2 is a pipeline for analysing data from sequenced wastewater samples and identifying given 
lineages of SARS-CoV-2. It was developed for wastewater samples, which are enriched for SARS-CoV-2. 
The pipeline can be used for continuous sampling. The output of the PiGx SARS-CoV-2 pipeline is summarized in a report
which provides an intuitive visual overview about the development of variant abundance over time and location 
(see our [example report](https://github.com/BIMSBbioinfo/pigx_sarscov2_ww#sample-reports)). 
Additionally there will be more detailed reports per sample, which cover the quality control of the samples, 
the detected variants and a taxonomic classification of all unaligned reads. This version of the pipeline was designed 
to work with paired-end amplicon sequencing data generated using the ARTIC nCoV-2019 primers.   

Features to enable e.g. single-end input are currently under development. We are happy to hear about any other feature 
request! Please consider using the Issue Tracker for this. 

# Installation

## Databases
Some tools require some databases to be installed locally. For this please either see the Documentation below or run 
[download_databases.sh](https://github.com/BIMSBbioinfo/pigx_sarscov2_ww/blob/main/scripts/download_databases.sh.in) after
the pipeline is installed. 

## Installation via Guix (recommended)

You can install this pipeline with all its dependencies using GNU Guix:
```sh 
guix install pigx-sars-cov2-ww
```

Using GNU Guix has many advantages in terms of reproducibility of your projects (see e.g. [here](https://academic.oup.com/gigascience/article/7/12/giy123/5114263)) 
If you don't have GNU Guix on your system you can also consider using [GNU Guix in a VM](https://guix.gnu.org/manual/en/html_node/Running-Guix-in-a-VM.html)


## Installation from source

You can also install PiGx-SARS-CoV-2 from source manually. Make sure that all the required dependencies (for this see e.g [configure.ac](https://github.com/BIMSBbioinfo/pigx_sarscov2_ww/blob/main/configure.ac)) 
are installed e.g. by installing them through a package manager like Conda. However, we can highly recommend using Guix (see above)

<details>
    <summary> The following tools must be available (latest version): </summary>

    - snakemake  
    - samtools  
    - bwa  
    - bedtools  
    - fastp  
    - fastqc  
    - R  
    - Rscript  
    - kraken2  
    - kraken2-build  
    - ktImportKrona  
    - ktImportTaxonomy  
    - ivar  
    - lofreq  
    - vep  
    - multiqc  
    - pandoc  

And the R-packages:  

    - DT  
    - base64url    
    - dplyr    
    - ggplot2    
    - magrittr  
    - plotly  
    - qpcR  
    - rmarkdown  
    - stringr  
    - tidyr  
    - reshape2  
    - R.utils  

All of these dependencies must be present in the environment at
configuration time.
</details>

Then you can use the buildsteps:  
```sh
./bootstrap.sh # to generate the "configure" script
./configure
make
make check
```

## Use the Docker Image
Will be available soon

# Documentation

Please see our documentation in order to find out about more details e.g. about the required structure of the input files:
https://bioinformatics.mdc-berlin.de/pigx_docs/pigx-sars-cov2-ww.html

# Sample Reports

See as example the HTML report produced by PiGx SARS-Cov-2 used for 
["COVID-19 infection dynamics revealed by SARS-CoV-2 wastewater sequencing analysis and deconvolution"](https://www.medrxiv.org/content/10.1101/2021.11.30.21266952v1)
https://bimsbstatic.mdc-berlin.de/akalin/AAkalin_pathogenomics/sarscov2_ww_reports/211104_pub_version/index.html 

# Getting help

If you have any questions please e-mail: pigx@googlegroups.com or use the web form to ask questions https://groups.google.com/forum/#!forum/pigx/. 

If you run into any bugs, please open an issue here: https://github.com/BIMSBbioinfo/pigx_sarscov2_ww/issues. 

# Links

- [The Bioinformatics and Omics Data Science Platform](https://bioinformatics.mdc-berlin.de)
