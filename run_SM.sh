numjobs=36


i=1
LOG="./nohup_SM_submit_"${i}".log"

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

