

<a name="logo"/>
<div align="center">
<img src="images/Logo_PIGx.png" alt="PiGx Logo"  width="30%" height="30%" ></img>
</a>
</div>

######  Copyright 2017: Alexander Gosdschan, Katarzyna Wreczycka, Bren Osberg.
######  This work is distributed under the terms of the GNU General Public License,  and is free to use for academic purposes

-----
# Summary

PIGx is a data processing pipeline for raw fastq read data of bisulfite experiments; it produces methylation and coverage information and can be used to produce information on differential methylation and segmentation. It was first developed by the Akalin group at MDC in Berlin in 2017.

The figure below provides a sketch of the process.
![](images/pipelineIO_BSseq.png )

# Getting started
PIGx can be installed through Guix (details to be added later) and consists of the following scripts: 

| File          | Purpose       |
| ------------- |:-------------:|
|  BSseq_pipeline.py | (main) Defines the rules for data processing.     |
| config.json   | Defines various parameters; e.g. input/output folder paths, sample names, etc. |
| func_defs.py  | Subscript that defines various functions called in the main snakemake script.                          |
| run_SM.sh     | Establishes indexed log files with time stamps and sets various parameters.      |

The dependencies are as follows:
 - FASTQC                        
 - TRIMGALORE                   
 - CUTADAPT                      
 - BISMARK_GENOME_PREPARATION    
 - BISMARK                       
 - BOWTIE2                       
 - DEDUPLICATE_BISMARK           
 - BISMARK_METHYLATION_EXTRACTOR 
 - BISMARK2REPORT                
 - SAMTOOLS 

All of these dependencies must be present in folder indicated in the config.json file by  ["paths"]["GTOOLBOX"].

There are additional dependencies that do need to be installed, but not located in ["paths"]["GTOOLBOX"]:

 - [python-rp2](https://rpy2.bitbucket.io/)
 - [pandoc](http://pandoc.org/)
 - [methylKit](https://github.com/al2na/methylKit)[>=1.3.1]
 - [genomation](http://bioinformatics.mdc-berlin.de/genomation/)
 - [DT](https://rstudio.github.io/DT/) 
 - [annotationhub](https://www.bioconductor.org/packages/release/bioc/html/AnnotationHub.html)
 - [rtracklayer](http://bioconductor.org/packages/release/bioc/html/rtracklayer.html)
 - [rmarkdown](http://rmarkdown.rstudio.com/)[>=1.5]

To run PIGx on your experimental data, first enter the necessary parameters in the spreadsheet file (see following section), and then from the terminal type

```
$ PIGx.x [options] 
```

---- 
# Input parameters

The input parameters specifying the desired behaviour of PIGx should be entered into the file Samplesheet.
When PIGx is run, the data from this file will be used to automatically generate a configuration file with the following values:
 
| Variable name | default value | description |
| ------------- |:-------------:|:-----------:|
| NICE          |    19         | integer: from -20 to 19; higher values make the program execution less demanding on computational resources |
| RCODE         |  ".read"      | string: describing the naming convention of files within the full data set  |
| directional   |    TRUE       | boolean: : is the sequencing experiment directional |
| NUMTHREADS    |     2         | integer: number of CPU threads the program should request for various functions that allow hyperthreading. |
|  INEXT        |  ".fq.gz"     | string:  input extension describes the standard suffix for filenames in a given data set  |
| SAMPLES       |     --        | struct: list of all the samples to be considered. |
|  files        |     --        | string(s): part of SAMPLE: lists files (without extension) to read. when 2 are specified, paired-end is assumed, otherwise, single end. |
|  PATHIN       |    "./"       | string: location of the experimental data files (.fastq[.gz|.bz2])   |
|  PATHOUT      |    "./out"    | string: ultimate location of the output data and report files   |
|  GENOMEPATH   |     ---       | string: location of the reference genome data for alignment   |
|  GTOOLBOX     | "~/.guix-profile/bin/"  | string: executable source files for the PIGx dependencies   |

 
All input files (paired or single end) must be present in the foler indicated by _PATHIN_, And must have their files listed among ["SAMPLES"]. Output from the snakemake script will then be sent to the folder indicated by ["paths"]["PATHOUT"], with subdirectories corresponding to the various stages of the process.

 
By default, the pipeline currently terminates at the stage of a final bismark report, however the 
main snakemake file still contains various comment lines that can be uncommented for the purpose 
of generating output from intermediate rules (i.e. without proceeding to this final stage)

The folder indicated by ["paths"]["GENOMEPATH"] must contain the reference genome being mapped to.

