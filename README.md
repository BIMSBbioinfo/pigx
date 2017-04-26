# pigx_bsseq: a pipeline for raw fastq read data of bisulfite experiments. | produces methylation and coverage information.
######  Developed by the Akalin group at MDC, Berlin, 2017

######  Copyright 2017: Alexander Gosdschan, Katarzyna Wreczycka, Bren Osberg.
######  This work is distributed under the terms of the GNU General Public License, 
######  and is free to use for academic purposes

-----

Various scripts are used in this package as listed in the table below:

| File          | Purpose       |
| ------------- |:-------------:|
|  BSseq_pipeline.py :|: Main script that defines the rules by which data conversion is performed in various steps      |
| config.json   | Config file that lists various parameters related to input/output folders, sample names, and programs |
| func_defs.py  | Subscript that defines various functions called in the main snakemake script                          |
| run_SM.sh     | Establishes indexed log files with time stamps and sets various parameters, such as the number of processes, etc.      |

- in addition to the data files themselves (generally compressed fastq)

----

The config file must specify the following features of the input data:

| Variable name | default value | description |
| NICE          |    19         | integer: from -20 to 19; higher values make the program execution less demanding on computational resources |
| RCODE         |  ".read"      | string: describing the naming convention of files within the full data set  |
| directional   |    TRUE       | boolean: : is the sequencing experiment directional |
| NUMTHREADS    |     2         | integer: number of CPU threads the program should request for various functions that allow hyperthreading. |
|  INEXT        |   ".fq.gz"    | string:  input extension describes the standard suffix for filenames in a given data set  |
| SAMPLES       |     --        | struct: list of all the samples to be considered. |
|  files        |     --        | string(s): part of SAMPLE: lists files (without extension) to read. when 2 are specified, paired-end is assumed, otherwise, single end. |
       
 
By default, the pipeline currently terminates at the stage of a final bismark report, however the 
main snakemake file still contains various comment lines that can be uncommented for the purpose 
of generating output from intermediate rules (i.e. without proceeding to this final stage)

Other commented lines may be useful for debugging, checking output, such as the interactive lines:
# import IPython;
# IPython.embed()

The snakemake script is dependent on the following functions:

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


all of which must be present in folder indicated in the config.json file by  ["paths"]["GTOOLBOX"]

All input files (paired or single end) must be present in the foler indicated by  ["paths"]["PATHIN"]
And must have their files listed in ["SAMPLES"]

Output from the snakemake script will then be sent to the folder indicated by ["paths"]["PATHOUT"], with 
subdirectories corresponding to the various stages of the process.

The folder indicated by ["paths"]["GENOMEPATH"] must contain the reference genome being mapped to.

