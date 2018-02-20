# ---------------------------------------------------------------------------- #
import glob
import fnmatch
import os
import re
import sys
import yaml

include: os.path.join(config['locations']['pkglibexecdir'], 'scripts/SnakeFunctions.py')
include: os.path.join(config['locations']['pkglibexecdir'], 'scripts/Check_Config.py')

localrules: makelinks

SCRIPT_PATH       = os.path.join(config['locations']['pkglibexecdir'], 'scripts/')
RULES_PATH        = os.path.join(config['locations']['pkglibexecdir'], 'Rules/')
SAMPLE_SHEET_FILE = config['locations']['sample-sheet']
REPORT_TEMPLATE   = os.path.join(SCRIPT_PATH,'Sample_Report.rmd')

# ---------------------------------------------------------------------------- #
# sample sheet input
with open(SAMPLE_SHEET_FILE, 'r') as stream:
    SAMPLE_SHEET = yaml.load(stream)

# ---------------------------------------------------------------------------- #
# check settings and sample_sheet validity
validate_config(config, SAMPLE_SHEET_FILE)

# ---------------------------------------------------------------------------- #
# Software executables
SOFTWARE = config['tools']


# Per sample software parameters:
# Loops throug the sample sheet and extracts per sample software parameters
# Flattens all samples with custom parameters into one dict - sample names must be unique
custom_param_names = sorted(list(set(SAMPLE_SHEET.keys()) - set(['samples','hub'])))
CUSTOM_PARAMS = dict()
for param_set in custom_param_names:
    for sample_name in SAMPLE_SHEET[param_set].keys():
        sample_set = SAMPLE_SHEET[param_set][sample_name]
        if isinstance(sample_set, dict) and 'params' in set(sample_set.keys()):
            CUSTOM_PARAMS[sample_name] = sample_set['params']

# ---------------------------------------------------------------------------- #
# Variable definition

# Default Function Parameters
PARAMS       = config['general']['params']
GENOME       = config['general']['assembly']
GENOME_ORIG  = config['locations']['genome-file']
PATH_FASTQ   = config['locations']['input-dir']
ANNOTATION   = config['locations']['gff-file']

# Sample name definition
PEAK_NAMES   = []
NAMES        = SAMPLE_SHEET['samples'].keys()

# Directory structure definition
OUTPUT_DIR      = config['locations']['output-dir']
# PATH_FASTQ      = os.path.join(OUTPUT_DIR, 'Fastq')
PATH_MAPPED     = os.path.join(OUTPUT_DIR, 'Mapped/Bowtie')
PATH_QC         = os.path.join(OUTPUT_DIR, 'FastQC')
PATH_INDEX      = os.path.join(OUTPUT_DIR, 'Bowtie2_Index')
PATH_LOG        = os.path.join(OUTPUT_DIR, 'Log')
PATH_PEAK       = os.path.join(OUTPUT_DIR, 'Peaks/MACS2')
PATH_BW         = os.path.join(OUTPUT_DIR, 'BigWig')
PATH_IDR        = os.path.join(OUTPUT_DIR, 'Peaks/IDR')
PATH_HUB        = os.path.join(OUTPUT_DIR, 'UCSC_HUB')
PATH_ANALYSIS   = os.path.join(OUTPUT_DIR, "Analysis")
PATH_ANNOTATION = os.path.join(OUTPUT_DIR, 'Annotation')
PATH_REPORTS    = os.path.join(OUTPUT_DIR, 'Reports')


# Directory structure for saved R objects
PATH_RDS            = os.path.join(PATH_ANALYSIS, 'RDS')
PATH_RDS_ANNOTATION = os.path.join(PATH_RDS, 'Annotation')
PATH_RDS_FEATURE    = os.path.join(PATH_RDS, 'Feature_Combination')
PATH_RDS_CHIPQC     = os.path.join(PATH_RDS, 'ChIPQC')
PATH_RDS_TEMP       = os.path.join(PATH_RDS, 'Temp')


# Hub variables which describe the types of files that can be used in the hub
TRACK_PATHS = {
    'bigWig' : {'path': PATH_MAPPED, 'suffix': 'bw', 'type':'bigWig'},
    'macs' :{'path': PATH_PEAK, 'suffix': 'bb', 'type':'bigBed'},
    'idr' : {'path': PATH_IDR, 'suffix': 'IDR.narrowPeak'}
}

# Collects the locations of all peaks
PEAK_NAME_LIST = {}


# names for markdown chunks which will be knit
# current names: ChIPQC, Extract_Signal_Annotation, Peak_Statistics, Annotate_Peaks
# ---------------------------------------------------------------------------- #


# Constructs the genome index prefix name
# If the prefix name is already supplied in the config file, it extracts this index
# instead of creating a new index
# index-dir should be the location of the index, without the genome prefix (i.e. without hg19)
if GENOME_ORIG == None:
    prefix_default = ''
    GENOME_ORIG = ''
else:
    prefix_default = os.path.join(PATH_INDEX, GENOME)

# set_default is a function of three arguments
# key, default, dictionary
# if dictionary[key] does not exist, it returns default
INDEX_PREFIX_NAME = set_default('index_prefix', GENOME,  config['general'])

# defines the link to the genome location
GENOME_PREFIX_PATH = os.path.join(prefix_default, INDEX_PREFIX_NAME)

# defines the location of the index (if the index is not specified)
INDEX_PREFIX_PATH  = os.path.join(set_default('index-dir', prefix_default, config['locations']), INDEX_PREFIX_NAME)



# ---------------------------------------------------------------------------- #
# due to the iditotic namig scheme in FASTQC the next lines construct
# FASTQC output files
FASTQC_DICT = {}
for i in NAMES:
    for fqfile in SAMPLE_SHEET['samples'][i]['fastq']:
        prefix  = fqfile
        prefix  = re.sub('.fq.*'   , '', prefix)
        prefix  = re.sub('.fastq.*', '', prefix)
        fastqc  = os.path.join(PATH_QC, i, prefix + "_fastqc.zip")
        FASTQC_DICT[prefix]  =  {'fastq'  : os.path.join(PATH_FASTQ, fqfile),
                                 'fastqc' : fastqc}
# ---------------------------------------------------------------------------- #
# RULE ALL
# Default output files from the pipeline

COMMAND         = []
GENOME_FASTA    = [GENOME_PREFIX_PATH + '.fa']
INDEX           = [INDEX_PREFIX_PATH  + '.1.bt2']
BOWTIE2         = expand(os.path.join(PATH_MAPPED, "{name}", "{name}.sorted.bam.bai"), name=NAMES)
BOWTIE2_STATS   = [os.path.join(PATH_RDS, "BowtieLog.rds")]
CHRLEN          = [GENOME_PREFIX_PATH + '.chrlen.txt']
TILLING_WINDOWS = [GENOME_PREFIX_PATH + '.GenomicWindows.GRanges.rds']
NUCLEOTIDE_FREQ = [GENOME_PREFIX_PATH + '.NucleotideFrequency.GRanges.rds']
FASTQC          = [FASTQC_DICT[i]['fastqc'] for i in list(FASTQC_DICT.keys())]
MULTIQC         = [os.path.join(PATH_REPORTS, "multiqc.html")]
ChIPQC          = expand(os.path.join(PATH_RDS_CHIPQC, "{name}_ChIPQC.rds"), name=NAMES)
BW              = expand(os.path.join(os.getcwd(), PATH_MAPPED, "{name}", "{name}.bw"),  name=NAMES)
LINKS           = expand(os.path.join(PATH_BW,  "{ex_name}.bw"),  ex_name=NAMES)


COMMAND = GENOME_FASTA + INDEX + BOWTIE2 + CHRLEN + TILLING_WINDOWS + NUCLEOTIDE_FREQ+ BW + LINKS + FASTQC + MULTIQC + ChIPQC + BOWTIE2_STATS


# ----------------------------------------------------------------------------- #
# include rules
include: os.path.join(RULES_PATH, 'Mapping.py')
include: os.path.join(RULES_PATH, 'Parse_Bowtie2log.py')
include: os.path.join(RULES_PATH, 'FastQC.py')
include: os.path.join(RULES_PATH, 'MultiQC.py')
include: os.path.join(RULES_PATH, 'ChIPQC.py')
include: os.path.join(RULES_PATH, 'BamToBigWig.py')


# ----------------------------------------------------------------------------- #
# ----------------------------------------------------------------------------- #
# CONDITIONAL OUTPUT FILES

peak_index = 'peak_calling' in set(SAMPLE_SHEET.keys())
if peak_index:
    if len(SAMPLE_SHEET['peak_calling'].keys()) > 0:
        PEAK_NAMES = SAMPLE_SHEET['peak_calling'].keys()
        MACS  = []
        QSORT = []
        suffix = 'narrowPeak'
        for name in PEAK_NAMES:
            suffix = get_macs2_suffix(name, CUSTOM_PARAMS)

            MACS    = MACS  + [os.path.join(PATH_PEAK,  name, name + "_peaks." + suffix)]
            QSORT   = QSORT + [os.path.join(PATH_PEAK,  name, name + "_qsort.bed" )]
            PEAK_NAME_LIST[name] = QSORT[-1]
        
        # ------------------------------------------------------------------------ #
        PEAK_STATISTICS = [os.path.join(PATH_RDS, "Peak_Statistics.rds")]

        include: os.path.join(RULES_PATH, 'Peak_Calling.py')
        include: os.path.join(RULES_PATH, 'Peak_Statistics.py')
        
        COMMAND = COMMAND + MACS + QSORT + PEAK_STATISTICS

# # ----------------------------------------------------------------------------- #
if 'idr' in set(SAMPLE_SHEET.keys()):
    if len(SAMPLE_SHEET['idr'].keys()) > 0:
        IDR = []
        for name in SAMPLE_SHEET['idr'].keys():
            IDR = IDR + [os.path.join(PATH_IDR, name, name + ".bed")]
            PEAK_NAME_LIST[name] = IDR[-1]

        include: os.path.join(RULES_PATH, 'IDR.py')
        COMMAND = COMMAND + IDR

# # ----------------------------------------------------------------------------- #
HUB_NAME = None
if 'hub' in set(SAMPLE_SHEET.keys()):
    HUB_NAME = SAMPLE_SHEET['hub']['name']
    HUB = [os.path.join(PATH_HUB, HUB_NAME, 'done.txt')]
    BB  = expand(os.path.join(PATH_PEAK,  "{name}", "{name}.bb"),  name=SAMPLE_SHEET['peak_calling'].keys())

    include: os.path.join(RULES_PATH, 'UCSC_Hub.py')
    COMMAND = COMMAND + BB + HUB


# ---------------------------------------------------------------------------- #
gtf_index = type(ANNOTATION) is str
if gtf_index:
    LINK_ANNOTATION    = [os.path.join(PATH_ANNOTATION, 'GTF_Link.gtf')]
    PREPARE_ANNOTATION = [os.path.join(PATH_ANNOTATION, 'Processed_Annotation.rds')]

    EXTRACT_SIGNAL_ANNOTATION = expand(os.path.join(PATH_RDS_TEMP,'{name}','{name}.Extract_Signal_Annotation.rds'), name=NAMES)
    
    # ------------------------------------------------------------------------ #
    # Formats the annotation + extracts signal profiles around pre-specified annotation regions
    include: os.path.join(RULES_PATH, 'Prepare_Annotation.py')
    include: os.path.join(RULES_PATH, 'Extract_Signal_Annotation.py')
    COMMAND = COMMAND + LINK_ANNOTATION + PREPARE_ANNOTATION + EXTRACT_SIGNAL_ANNOTATION

    if peak_index:
        ANNOTATE_PEAKS     = expand(os.path.join(PATH_RDS_TEMP,'{name}','{name}.Annotate_Peaks.rds'), name=PEAK_NAMES)
        include: os.path.join(RULES_PATH, 'Annotate_Peaks.py')
        COMMAND = COMMAND + LINK_ANNOTATION + PREPARE_ANNOTATION + ANNOTATE_PEAKS

# ---------------------------------------------------------------------------- #
if 'feature_combination' in set(SAMPLE_SHEET.keys()):
    FEATURE_NAMES = SAMPLE_SHEET['feature_combination'].keys()
    if len(FEATURE_NAMES) > 0:

        FEATURE = expand(os.path.join(PATH_RDS_FEATURE,'{name}_FeatureCombination.rds'),
        name = FEATURE_NAMES)

    include: os.path.join(RULES_PATH, 'Feature_Combination.py')
    COMMAND = COMMAND + FEATURE

# ----------------------------------------------------------------------------- #
# REPORT INPUT
REPORT         = [os.path.join(PATH_REPORTS, 'ChIP_Seq_Report.html')]
# REPORT_CHUNKS  = ['EXTRACT_SIGNAL_ANNOTATION','PEAK_STATISTICS','ANNOTATE_PEAKS','ChIPQC']
# This lines convert the analysis code chunks from SNAKEMAKE rule language to R language
REPORT_CHUNKS  = {'EXTRACT_SIGNAL_ANNOTATION':'Extract_Signal_Annotation','PEAK_STATISTICS':'Peak_Statistics','ANNOTATE_PEAKS':'Annotate_Peaks','ChIPQC':'ChIPQC'}
REPORT_INPUT   = []
ANALISYS_NAMES = []
for i in REPORT_CHUNKS.keys():
    if i in globals().keys():
        REPORT_INPUT   = REPORT_INPUT   + globals()[i]
        ANALISYS_NAMES = ANALISYS_NAMES + [i]

ANALISYS_NAMES = [REPORT_CHUNKS[i] for i in ANALISYS_NAMES]
include: os.path.join(RULES_PATH, 'Knit_Report.py')
COMMAND = COMMAND + REPORT
# ----------------------------------------------------------------------------- #
rule all:
    input:
        COMMAND
