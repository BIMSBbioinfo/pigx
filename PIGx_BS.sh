#!/bin/bash
if [ $# -eq 0 ]
  then
    echo "No argument supplied. Please specify whether you wish to perform a "dry" run, a "cluster" submission, an "unlock"-ing of the snakemake directory, or a "local" calculation on the current machine. Exiting."
    exit
fi


#===== PATHS ===== #
runparams="scripts/PIGx_runtime_params.csv"
path2configfile="./config.json"
path2programsJSON="test_dataset/PROGS.json"

#========================================================================================
#----------  CREATE CONFIG FILE:  ----------------------------------------------

scripts/create_configfile.py $runparams $path2configfile $path2programsJSON


path_IN=$(         grep  PATHIN     PIGx_BS.input  | sed 's/.*"\(.*\)"/\1/g')
path_OUT=$(        grep  PATHOUT    PIGx_BS.input  | sed 's/.*"\(.*\)"/\1/g')
path_refGenomeIN=$(grep GENOMEPATH  PIGx_BS.input  | sed 's/.*"\(.*\)"/\1/g')

mkdir -p ${path_OUT}
mkdir -p ${path_OUT}"/path_links"


# create symbolic link in output directory linking directly to the source data
ln -sfn ${path_IN} ${path_OUT}"/path_links/input"
# likewise for the genome reference:
ln -sfn ${path_refGenomeIN} ${path_OUT}"/path_links/refGenome"
# N.B. Any previous links of the same name are over-written.


# ================================================================
# create a unique log file for this run:
i=1
LOG="./PIGx_submission_"${i}".log"
while [ -f ${LOG} ]
do
  i=$((i+1))
  LOG="./PIGx_submission_"${i}".log"
done
# ================================================================

echo "starting Snakemake session on " $(date) >  ${LOG}
echo "" >>${LOG}

echo ""                                       >> ${LOG}
echo "------ using samples : ------"          >> ${LOG}
grep -i "files" config.json                   >> ${LOG}

echo ""                                       >> ${LOG}
echo "------ from folder : --------"          >> ${LOG}
grep   "PATHIN" config.json                   >> ${LOG}

echo ""                                       >> ${LOG}
echo "------ to folder : ----------"          >> ${LOG}
grep   "PATHOUT" config.json                  >> ${LOG}

echo ""                                       >> ${LOG}

#========================================================================================
#----------  NOW START RUNNING SNAKEMAKE:  ----------------------------------------------
#TODO: do we need the "--forceall" option?

echo "------ Commencing Snakemake : ------"   >> ${LOG}

if [ $1 == "dry" ] 
then
snakemake -s BSseq_pipeline.py   --configfile $path2configfile -d $path_OUT  -n
elif [ $1 == "cluster" ] 
then
snakemake -s BSseq_pipeline.py   --configfile $path2configfile -d $path_OUT --jobs ${numjobs}  --cluster "qsub -V -l h_vmem=${MEM} -pe smp ${NUMTHREADS}  -l h_rt=${RUNTIME} -l h_stack=128k"   --latency-wait ${LATWAIT}  >> ${LOG} 2>&1  &
elif [ $1 == "local" ] 
then
snakemake -s BSseq_pipeline.py   --configfile $path2configfile -d $path_OUT  >> ${LOG} 
elif [ $1 == "unlock" ] 
then
snakemake -s BSseq_pipeline.py --unlock
else
echo "command line arg not understood. Exiting without submission"
fi
