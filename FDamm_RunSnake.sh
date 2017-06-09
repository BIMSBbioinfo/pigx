# ---------------------------------------------------------------------------------- #
SNAKEFILE='/home/vfranke/Projects/FDamm_scRNAseq/Scripts/Snakemake/Snake_Dropseq.py'
WORKDIR='/data/akalin/Projects/FDamm_scRNAseq'
snakemake -R --snakefile $SNAKEFILE --directory $WORKDIR --timestamp --jobs 24 --cluster "qsub -V -l h_vmem={params.mem} -pe smp {params.threads} -l h_rt=36:00:00" --rerun-incomplete --latency-wait 30
