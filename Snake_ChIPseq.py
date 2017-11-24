"""
SNAKEFILE='/home/vfranke/Projects/AAkalin_PIX/Snake_ChIPseq.py'
WORKDIR='/data/akalin/vfranke/AAkalin_PIX/ChIP'
CONFIGFILE='/home/vfranke/Projects/AAkalin_PIX/config.yaml'
PATH='/home/vfranke/bin/Software/miniconda3/envs/p35/bin:/usr/local/bin:/usr/bin:/bin:/home/vfranke/.guix-profile/bin:/home/vfranke/.guix-profile/sbin:/home/vfranke/bin'
# Beast run

snakemake -R --snakefile $SNAKEFILE --directory $WORKDIR --jobs 4 --rerun-incomplete --configfile $CONFIGFILE --latency-wait 30 --dryrun

snakemake -R --snakefile $SNAKEFILE --directory $WORKDIR --jobs 4 --rerun-incomplete --configfile $CONFIGFILE --latency-wait 30

snakemake -R --snakefile $SNAKEFILE --directory $WORKDIR --jobs 4 --rerun-incomplete --configfile $CONFIGFILE --latency-wait 30 --dryrun

# Visualize the graph horizontally
snakemake --dag --snakefile $SNAKEFILE --directory $WORKDIR --configfile $CONFIGFILE | perl -pe 's|graph\[|graph[rankdir=LR,|'| dot -T png > Snake_ChIPseq.png

# Cluster run
snakemake -R --snakefile $SNAKEFILE --directory $WORKDIR --timestamp --jobs 24 --cluster "qsub -V -l h_vmem={params.mem} -pe smp {params.threads} -l h_rt=36:00:00" --rerun-incomplete --latency-wait 30

for i in `ls $WORKDIR | grep -v Fastq$ | grep -v Genome`;do rm -r $WORKDIR/$i;done

rm $(snakemake --snakefile $SNAKEFILE --directory $WORKDIR --configfile $CONFIGFILE --summary | tail -n+2 | cut -f1)


gp -i gcc-toolchain@5 gfortran@5 pkg-config curl cairo libxt libxml2 openssl xz zlib -p pix

gp -i gcc -p pix

bash
cd /home/vfranke/bin/guix
profile=/home/vfranke/bin/guix/pix
source $profile/etc/profile

export PATH="$profile/bin${PATH:+:}$PATH"
export PKG_CONFIG_PATH="pix/lib/pkgconfig${PKG_CONFIG_PATH:+:}$PKG_CONFIG_PATH"
export XDG_DATA_DIRS="pix/share${XDG_DATA_DIRS:+:}$XDG_DATA_DIRS"
export GIO_EXTRA_MODULES="pix/lib/gio/modules${GIO_EXTRA_MODULES:+:}$GIO_EXTRA_MODULES"
export R_LIBS_SITE="pix/site-library/${R_LIBS_SITE:+:}$R_LIBS_SITE"
export GUIX_LOCPATH=$profile/lib/locale
export CURL_CA_BUNDLE="/etc/ssl/certs/ca-bundle.crt"
/home/vfranke/bin/guix/pix/bin/R
source("https://bioconductor.org/biocLite.R")
biocLite('S4Vectors')
biocLite(c('data.table','genomation','GenomicRanges','RCAS','rtracklayr','stringr', 'plotly', 'DT', 'data.table', 'topGO', 'motifRG', 'biomaRt', 'AnnotationDbi', 'GenomicRanges', 'BSgenome.Hsapiens.UCSC.hg19', 'GenomeInfoDb', 'Biostrings', 'rtracklayer', 'org.Hs.eg.db', 'GenomicFeatures', 'genomation', 'rmarkdown', 'knitr', 'S4Vectors'), dependencies=TRUE)




Integrate_Expression_Mouse_Hamster_Difference
1. Figure out how to prevent bombing

### 13. config file check

Q:
1. How to write dynamic rules? - that the rules execute based on the input parameters?
    - if the user supplies the genome index, do not make the index
2. What to do downstream? Some basic plots? ChromHMM
3. How to speed up the extension?

"""

SCRIPT_PATH = os.path.join(workflow.basedir,'Scripts')
RULES_PATH  = os.path.join(workflow.basedir,'Rules')
PARAMS_PATH = '.'

# ---------------------------------------------------------------------------- #
import glob
import fnmatch
import os
import re
import sys
import yaml

# from SnakeFunctions import *
include: 'SnakeFunctions.py'
from Check_Config import *
localrules: makelinks


# ---------------------------------------------------------------------------- #
# check config validity
if check_config(config) == 1:
    quit()

# ---------------------------------------------------------------------------- #
# Software parameters
with open(os.path.join(workflow.basedir, 'software.yaml'), 'r') as stream:
    SOFTWARE_CONFIG = yaml.load(stream)

# Function parameter
APP_PARAMS = SOFTWARE_CONFIG['params']

# Software executables
SOFTWARE   = SOFTWARE_CONFIG['software']
# ---------------------------------------------------------------------------- #
# Variable definition
GENOME       = config['genome']['name']
GENOME_FASTA = config['genome']['fasta']
NAMES        = config['samples'].keys()
PATH_FASTQ   = config['fastq']
PARAMS       = config['params']
ANNOTATION   = config['annotation']


# Directory structure definition
PATH_MAPPED     = 'Mapped/Bowtie'
PATH_QC         = 'FastQC'
PATH_INDEX      = 'Bowtie2_Index'
PATH_LOG        = 'Log'
PATH_PEAK       = 'Peaks/MACS2'
PATH_BW         = 'BigWig'
PATH_IDR        = 'IDR'
PATH_HUB        = 'UCSC_HUB'
PATH_ANALYSIS   = "Analysis"

# Analisys path
PATH_FEATURE    = os.path.join(PATH_ANALYSIS, 'Feature_Combination')

# Directory structure for saved R objects
PATH_RDS            = os.path.join(PATH_ANALYSIS, 'RDS')
PATH_RDS_ANNOTATION = os.path.join(PATH_RDS, 'Annotation')
PATH_RDS_ANALYSIS   = os.path.join(PATH_RDS, 'Analysis')


# Hub variables which describe the types of files that can be used in the hub
TRACK_PATHS = {
    'bigWig' : {'path': PATH_MAPPED, 'suffix': 'bw', 'type':'bigWig'},
    'macs' :{'path': PATH_PEAK, 'suffix': 'bb', 'type':'bigBed'},
    'idr' : {'path': PATH_IDR, 'suffix': 'IDR.narrowPeak'}
}


#
# ---------------------------------------------------------------------------- #
# config defaults
if not ('extend' in config['params'].keys()):
    config['params']['extend'] = 0


if GENOME_FASTA == None:
    prefix_default = ''
    GENOME_FASTA = ''
else:
    prefix_default = os.path.join(PATH_INDEX, GENOME)


PREFIX = os.path.join(set_default('index', prefix_default, config), GENOME)
print(PREFIX)

# ----------------------------------------------------------------------------- #
# include rules
include: os.path.join(RULES_PATH, 'Mapping.py')
include: os.path.join(RULES_PATH, 'FastQC.py')
include: os.path.join(RULES_PATH, 'BamToBigWig.py')

# ----------------------------------------------------------------------------- #
# RULE ALL
COMMAND    = []
INDEX      = [PREFIX + '.1.bt2']
BOWTIE2    = expand(os.path.join(PATH_MAPPED, "{name}", "{name}.sorted.bam.bai"), name=NAMES)
CHRLEN     = [PREFIX + '.chrlen.txt']
FASTQC     = expand(os.path.join(PATH_QC,     "{name}", "{name}.fastqc.done"), name=NAMES)
BW         = expand(os.path.join(os.getcwd(), PATH_MAPPED, "{name}", "{name}.bw"),  name=NAMES)
LINKS      = expand(os.path.join(PATH_BW,  "{ex_name}.bw"),  ex_name=NAMES)

COMMAND = COMMAND + INDEX + BOWTIE2 + CHRLEN + BW + LINKS + FASTQC



# ----------------------------------------------------------------------------- #
if 'peak_calling' in set(config.keys()):
    MACS  = []
    QSORT = []
    suffix = 'narrowPeak'
    for name in config['peak_calling'].keys():
        suffix = get_macs2_suffix(name, config)

        MACS    = MACS  + [os.path.join(PATH_PEAK,  name, name + "_peaks." + suffix)]
        QSORT   = QSORT + [os.path.join(PATH_PEAK,  name, name + "_qsort." + suffix)]

    include: os.path.join(RULES_PATH, 'Peak_Calling.py')
    COMMAND = COMMAND + MACS + QSORT

# # ----------------------------------------------------------------------------- #
# if 'idr' in set(config.keys()):
    IDR     = expand(os.path.join(PATH_IDR,    "{name}", "{name}.narrowPeak"),     name=config['idr'].keys())

    include: os.path.join(RULES_PATH, 'IDR.py')
    COMMAND = COMMAND + IDR

# # ----------------------------------------------------------------------------- #
HUB_NAME = None
if 'hub' in set(config.keys()):
    HUB_NAME = config['hub']['name']
    HUB = [os.path.join(PATH_HUB, HUB_NAME, 'done.txt')]
    BB  = expand(os.path.join(PATH_PEAK,  "{name}", "{name}.bb"),  name=config['peak_calling'].keys())

    include: os.path.join(RULES_PATH, 'UCSC_Hub.py')
    COMMAND = COMMAND + BB + HUB


# # ----------------------------------------------------------------------------- #
if 'feature_combination' in set(config.keys()):
    peak_files = []
    if 'idr' in config['feature_combination'].keys():
        peak_files = peak_files + [os.path.join(PATH_IDR,  x, x + '.narrowPeak') for x in config['feature_combination']['idr']]

    if 'peaks' in config['feature_combination'].keys():
        peak_files = peak_files + [os.path.join(PATH_PEAK, x, x + '_peaks.narrowPeak') for x in config['feature_combination']['peaks']]
    FEATURE = [os.path.join(PATH_FEATURE,'Feature_Combination.tsv')]


    include: os.path.join(RULES_PATH, 'Feature_Combination.py')
    COMMAND = COMMAND + FEATURE
#
#
# # ----------------------------------------------------------------------------- #
if 'gtf' in set(config['annotation'].keys()):
    PREPARE_ANNOTATION = [os.path.join(PATH_RDS_ANNOTATION,'Processed_Annotation.rds')]

    ANNOTATE_PEAKS     = expand(os.path.join(PATH_RDS_ANALYSIS,'{name}','{name}.Annotate_Peaks.rds'), name=NAMES)

    EXTRACT_SIGNAL_PEAKS = expand(os.path.join(PATH_RDS_ANALYSIS,'{name}','{name}.Extract_Signal_Peaks.rds'), name=NAMES)

    EXTRACT_SIGNAL_ANNOTATION = expand(os.path.join(PATH_RDS_ANALYSIS,'{name}','{name}.Extract_Signal_Annotation.rds'), name=NAMES)


    # include: os.path.join(RULES_PATH, 'Prepare_Annotation.py')
    # include: os.path.join(RULES_PATH, 'Extract_Signal_Annotation.py')
    # COMMAND = COMMAND + PREPARE_ANNOTATION + ANNOTATE_PEAKS + EXTRACT_SIGNAL_PEAKS + EXTRACT_SIGNAL_ANNOTATION

# ----------------------------------------------------------------------------- #
# if 'feature_combination' in set(config.keys()) and 'gtf' in set(config['annotation'].keys()





# ----------------------------------------------------------------------------- #
rule all:
    input:
        COMMAND
