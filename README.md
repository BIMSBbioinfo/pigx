# SNAKEMAKE PIPELINE DEVELOPED BY WRITTEN BY THE AKALIN GROUP AT MDC, BERLIN, 2017
# ---last updated on  Wed Mar 22 13:12:14 CET 2017  by  blosberg  at location  , Bren-Osbergs-MacBook.local

#  changes from  Wed Mar 22 13:12:14 CET 2017 : first README text description provided.


# For processing bisulfite sequencing data from raw fastq files 

# Copyright 2017: Alexander Gosdschan, Katarzyna Wreczycka, Bren Osberg.
# This work is distributed under the terms of the GNU General Public License, 
# and is free to use for academic purposes
#================================================================================================

The main snakemake executable is BSseq_pipeline.py; it requires the following input files:

* config.json   #---- a config file that lists various parameters related to input/output folders, sample names, and programs
* func_defs.py  #---- a subscript that defines various functions called in the main snakemake script.

in addition to the data files themselves (generally compressed fastq)

The script itself can be called from a bash executable 

* run_SM.sh     #---- which establishes indexed log files with time stamps and sets various parameters, such as the number of processes, etc.


By default, the pipeline currently terminates at the stage of a final bismark report, however the 
main snakemake file still contains various comment lines that can be uncommented for the purpose 
of generating output from intermediate rules (i.e. without proceeding to this final stage

Other commented lines may be useful for debugging, checking output, such as the interactive lines:

# import IPython;
# IPython.embed()

The snake make script is dependent on the following functions:

FASTQC                        
TRIMGALORE                   
CUTADAPT                      
BISMARK_GENOME_PREPARATION    
BISMARK                       
BOWTIE2                       
DEDUPLICATE_BISMARK           
BISMARK_METHYLATION_EXTRACTOR 
BISMARK2REPORT                
SAMTOOLS 


all of which must be present in folder indicated in the config.json file by  ["paths"]["GTOOLBOX"]

All input files (paired or single end) must be present in the foler indicated by  ["paths"]["PATHIN"]
And must have their files listed in ["SAMPLES"]

Output from the snakemake script will then be sent to the folder indicated by ["paths"]["PATHOUT"], with 
subdirectories corresponding to the various stages of the process.

The folder indicated by ["paths"]["GENOMEPATH"] must contain the reference genome being mapped to.

                    
