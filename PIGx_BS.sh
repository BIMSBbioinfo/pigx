#!/bin/bash

#===== DEFAULT PATHS ===== #
tablesheet="test_dataset/TableSheet_test.csv"
path2configfile="./config.json"
path2programsJSON="test_dataset/PROGS.json"

#=========== PARSE PARAMETERS ============#

usage="

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


DESCRIPTION


PIGx is a data processing pipeline for raw fastq read data of bisulfite experiments.
It produces methylation and coverage information and can be used to produce information 
on differential methylation and segmentation. 

It was first developed by the Akalin group at MDC in Berlin in 2017.


USAGE: $(basename "$0") [-h] [-t|--tablesheet FILE] 
                        [-p|--programs FILE] 
                        [-c|--configfile FILE] 
                        [-s|--snakeparams PARAMS]

OPTIONAL ARGUMENTS: 

[-t|--tablesheet FILE]      The tablesheet containing the basic configuration information 
                            for running the BSseq_pipeline. 
                        
[-p|--programs FILE]        A json file containing the paths to the required tools.     

[-c|--configfile FILE]      The config file used for calling the underlying snakemake process.
                            By default the file '${path2configfile}' is dynamically created from tablesheet and 
                            programs file.

[-s|--snakeparams PARAMS]   Additional parameters to be passed down to snakemake, e.g. 
                                --dryrun    do not exectute anything
                                --forceall  rerun the whole pipeline 


"

# https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash
#
# Use -gt 1 to consume two arguments per pass in the loop (e.g. each
# argument has a corresponding value to go with it).
# Use -gt 0 to consume one or more arguments per pass in the loop (e.g.
# some arguments don't have a corresponding value to go with it such
# as in the --default example).
# note: if this is set to -gt 0 the /etc/hosts part is not recognized ( may be a bug )
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -t|--tablesheet)
    tablesheet="$2"
    shift # past argument
    ;;
    -c|--configfile)
    path2configfile="$2"
    shift # past argument
    ;;
    -p|--programs)
    path2programsJSON="$2"
    shift # past argument
    ;;
    -s|--snakeparams)
    snakeparams="$2"
    shift # past argument=value
    ;;
    -h|--help)
    echo "$usage"
    shift # past argument=value
    exit 1
    ;;
    # --default)
    # DEFAULT=YES
    # ;;
    *)
            # unknown option
    ;;
esac
shift # past argument or value
done

# echo "${tablesheet} ${path2configfile} ${path2programsJSON} ${snakeparams}" 

#======================================================================================
#----------  CREATE CONFIG FILE:  ----------------------------------------------

scripts/create_configfile.py $tablesheet $path2configfile $path2programsJSON

#======================================================================================
#----------  NOW CREATE SYMBOLIC LINKS TO THE INPUTS AND REFERENCE GENOME -------------

path_OUT=$( python -c "import sys, json; print(json.load(sys.stdin)['PATHOUT'])" < $path2configfile)
path_IN=$( python -c "import sys, json; print(json.load(sys.stdin)['PATHIN'])" < $path2configfile)
path_refG=$( python -c "import sys, json; print(json.load(sys.stdin)['GENOMEPATH'])" < $path2configfile)

mkdir -p ${path_OUT}
mkdir -p ${path_OUT}"path_links"
mkdir -p ${path_OUT}"path_links/input"

# create links within the output folder that point directly to the 
# reference genome, as well as to each sample input file  
# so that it's clear where the source data came from.
# N.B. Any previous links of the same name are over-written.

# link to reference genome:
ln -sfn ${path_refG} ${path_OUT}"/path_links/refGenome"

# create file links:
scripts/create_file_links.py $path2configfile 


#========================================================================================
#----------  NOW START RUNNING SNAKEMAKE:  ----------------------------------------------

pathout=$( python -c "import sys, json; print(json.load(sys.stdin)['PATHOUT'])" < $path2configfile)

snakemake -s BSseq_pipeline.py --configfile $path2configfile -d $pathout
