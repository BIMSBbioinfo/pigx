# ---------------------------------------------------------------------------------- #
SNAKEFILE='/home/vfranke/Projects/CBirchmeier_Neuro/Scripts/Snake/Snake_ChIPseq.py'
WORKDIR='/data/akalin/Projects/CBirchmeier_Neuro/Data/AccessoryData/Analysis'
snakemake -R --snakefile $SNAKEFILE --directory $WORKDIR --timestamp --jobs 24 --cluster "qsub -V -l h_vmem={params.mem} -pe smp {params.threads} -l h_rt=36:00:00" --rerun-incomplete --latency-wait 30 

# ---------------------------------------------------------------------------------- #
# in house data
SNAKEFILE='/home/vfranke/Projects/CBirchmeier_Neuro/Scripts/Snake/Snake_InHouse.py'
WORKDIR='/data/akalin/Projects/CBirchmeier_Neuro/Data/mm10_Insm1-Ascl1'
snakemake -R --snakefile $SNAKEFILE --directory $WORKDIR --timestamp --jobs 24 --cluster "qsub -V -l h_vmem={params.mem} -pe smp {params.threads} -l h_rt=36:00:00" --rerun-incomplete --latency-wait 30
