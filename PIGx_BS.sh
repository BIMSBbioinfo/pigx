#!/bin/bash

#===== PATHS ===== #
tablesheet="test_dataset/TableSheet_test.csv"
path2configfile="./config.json"
path2programsJSON="test_dataset/PROGS.json"

#========================================================================================
#----------  CREATE CONFIG FILE:  ----------------------------------------------

scripts/create_configfile.py $tablesheet $path2configfile $path2programsJSON


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

echo "------ Commencing Snakemake : ------"   >> ${LOG}

snakemake -s BSseq_pipeline.py --configfile $path2configfile -d $pathout
