#!/usr/bin/bash

additional_args=$1

echo "ARGS: ${additional_args}"

SRCDIR="/data/akalin/buyar/collaborations/jonathan/pipeline/pigx_crispr"

#echo "Removing previously generated output folder"
#rm -rf ${SRCDIR}/sample_data/output

settings="${SRCDIR}/sample_data/settings.yaml"
snakefile="${SRCDIR}/pigx_crispr.py"
snakemake='/home/buyar/.guix-profile/bin/snakemake'


${snakemake} ${additional_args} -p --configfile ${settings}  --snakefile ${snakefile} -j 4 
#${snakemake} --forceall -p --configfile ${settings}  --snakefile ${snakefile} -j 4 --dryrun
#	${snakemake} --forceall -p --configfile ${settings}  --snakefile ${snakefile} -j 4
