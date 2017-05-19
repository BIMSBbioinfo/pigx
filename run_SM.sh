numjobs=1

i=1
LOG="./nohup_SM_submit_"${i}".log"
# ================================================================

while [ -f ${LOG} ]
do
  i=$((i+1))
  LOG="./nohup_SM_submit_"${i}".log"
done

echo "starting Snakemake session on " $(date) >  ${LOG}
echo "" >>${LOG}

echo ""                                       >> ${LOG}
echo "------ using samples : ------"          >>${LOG}
grep -i "files" config.json                   >> ${LOG}

echo ""                                       >> ${LOG}
echo "------ from folder : --------"          >>${LOG}
grep   "PATHIN" config.json                   >> ${LOG}

echo ""                                       >> ${LOG}
echo "------ to folder : ----------"          >>${LOG}
grep PATHOUT config.json                      >> ${LOG}

echo ""                                       >> ${LOG}
echo "------ Commencing Snakemake : ------"   >>${LOG}

nohup  snakemake -s BSseq_pipeline.py --jobs ${numjobs}        >> ${LOG} &

# The following line can be used for cluster submission of the SM script:
# snakemake -s BSseq_pipeline_v2.py --jobs ${numjobs}   --cluster "qsub -V -l h_vmem=500M -pe smp 1  -l h_rt=32:00:00 -l h_stack=128k"       >> ${LOG} &
