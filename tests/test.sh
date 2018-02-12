echo "start: `date`"

srcdir='/data/local/buyar/pigx/scrna/pigx_scrnaseq/'

SNAKEFILE="${srcdir}/Snake_Dropseq.py"
WORKDIR="${srcdir}/tests/output"
CONFIGFILE="${srcdir}/tests/test.config.yaml"
PATH='/home/buyar/.conda/envs/pigx_scrna/bin:/usr/bin:/home/buyar/.guix-profile/bin/'

export PATH=$PATH
export TMPDIR='/data/local/buyar/tmp'
export _JAVA_OPTIONS=-Djava.io.tmpdir=$TMPDIR

snakemake -R -p --snakefile $SNAKEFILE --directory $WORKDIR --jobs 12 --rerun-incomplete --configfile $CONFIGFILE --latency-wait 30 --dryrun
snakemake -R -p --snakefile $SNAKEFILE --directory $WORKDIR --jobs 11 --rerun-incomplete --configfile $CONFIGFILE --latency-wait 30 

if ! test -f "${srcdir}/tests/output/Mapped/hg19.scRNA-Seq.report.html"
then
  echo "ERROR: could not find report output"
  exit 1
fi

echo "end: `date`"
