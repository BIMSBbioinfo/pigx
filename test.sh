dry=$1
SRCDIR="/data/akalin/buyar/collaborations/jonathan/pipeline/pigx_crispr"
settings="${SRCDIR}/sample_data/settings.yaml"
snakefile="${SRCDIR}/pigx_crispr.py"
snakemake='/home/buyar/.guix-profile/bin/snakemake'

if [ ${dry} == 'test' ] 
then
	${snakemake} -p --configfile ${settings}  --snakefile ${snakefile} -j 4 --dryrun
else
	${snakemake} -p --configfile ${settings}  --snakefile ${snakefile} -j 4
fi
