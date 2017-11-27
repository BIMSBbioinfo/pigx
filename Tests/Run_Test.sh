SNAKEFILE='/home/vfranke/Projects/AAkalin_PIX/Snake_ChIPseq.py'
WORKDIR='/data/akalin/vfranke/AAkalin_PIX/ChIP'
CONFIGFILE='/home/vfranke/Projects/AAkalin_PIX/config.yaml'
PATH='/home/vfranke/bin/Software/miniconda3/envs/p35/bin:/usr/local/bin:/usr/bin:/bin:/home/vfranke/.guix-profile/bin:/home/vfranke/.guix-profile/sbin:/home/vfranke/bin'
# Beast run

snakemake -R --snakefile $SNAKEFILE --directory $WORKDIR --jobs 4 --rerun-incomplete --configfile $CONFIGFILE --latency-wait 30 --dryrun

# snakemake -R --snakefile $SNAKEFILE --directory $WORKDIR --jobs 4 --rerun-incomplete --configfile $CONFIGFILE --latency-wait 30