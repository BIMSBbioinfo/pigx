# ---------------------------------------------------------------------------- #
import glob
import fnmatch
import os
import re
import sys
import yaml

# from SnakeFunctions import *
include: 'SnakeFunctions.py'
include: os.path.join(config['locations']['pkglibexecdir'], 'scripts/Check_Conf.py')

localrules: makelinks

SCRIPT_PATH       = os.path.join(config['locations']['pkglibexecdir'], 'scripts/')
RULES_PATH        = os.path.join(config['locations']['pkglibexecdir'], 'Rules/')
SAMPLE_SHEET_FILE = config['locations']['sample-sheet']

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
GENOME_FASTA = config['locations']['genome-file']
PATH_FASTQ   = config['locations']['input-dir']
ANNOTATION   = config['locations']['annotation']

# Sample name definition
PEAK_NAMES   = SAMPLE_SHEET['peak_calling'].keys()
NAMES        = SAMPLE_SHEET['samples'].keys()


# Directory structure definition
PATH_MAPPED     = 'Mapped/Bowtie'
PATH_QC         = 'FastQC'
PATH_INDEX      = 'Bowtie2_Index'
PATH_LOG        = 'Log'
PATH_PEAK       = 'Peaks/MACS2'
PATH_BW         = 'BigWig'
PATH_IDR        = 'Peaks/IDR'
PATH_HUB        = 'UCSC_HUB'
PATH_ANALYSIS   = "Analysis"
PATH_ANNOTATION = 'Annotation'


# Directory structure for saved R objects
PATH_RDS            = os.path.join(PATH_ANALYSIS, 'RDS')
PATH_RDS_ANNOTATION = os.path.join(PATH_RDS, 'Annotation')
PATH_RDS_FEATURE    = os.path.join(PATH_RDS, 'Feature_Combination')
PATH_RDS_TEMP       = os.path.join(PATH_RDS, 'Temp')


# Hub variables which describe the types of files that can be used in the hub
TRACK_PATHS = {
    'bigWig' : {'path': PATH_MAPPED, 'suffix': 'bw', 'type':'bigWig'},
    'macs' :{'path': PATH_PEAK, 'suffix': 'bb', 'type':'bigBed'},
    'idr' : {'path': PATH_IDR, 'suffix': 'IDR.narrowPeak'}
}

# Collects the locations of all peaks
PEAK_NAME_LIST = {}

#
# ---------------------------------------------------------------------------- #
# config defaults
if not ('extend' in PARAMS.keys()):
    PARAMS['extend'] = 0


# Constructs the genome index prefix name
# If the prefix name is already supplied in the config file, it extracts this index
# instead of creating a new index
# index-dir should be the location of the index, without the genome prefix (i.e. without hg19)
if GENOME_FASTA == None:
    prefix_default = ''
    GENOME_FASTA = ''
else:
    prefix_default = os.path.join(PATH_INDEX, GENOME)

INDEX_PREFIX_NAME = set_default('index_prefix', GENOME,  config['general'])
PREFIX = os.path.join(set_default('index-dir', prefix_default, config['locations']), INDEX_PREFIX_NAME)
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
if 'peak_calling' in set(SAMPLE_SHEET.keys()):
    if len(SAMPLE_SHEET['peak_calling'].keys()) > 0:
        MACS  = []
        QSORT = []
        suffix = 'narrowPeak'
        for name in PEAK_NAMES:
            suffix = get_macs2_suffix(name, CUSTOM_PARAMS)

            MACS    = MACS  + [os.path.join(PATH_PEAK,  name, name + "_peaks." + suffix)]
            QSORT   = QSORT + [os.path.join(PATH_PEAK,  name, name + "_qsort.bed" )]
            PEAK_NAME_LIST[name] = QSORT[-1]
        
        include: os.path.join(RULES_PATH, 'Peak_Calling.py')
        COMMAND = COMMAND + MACS + QSORT

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
gtf_index = 'gtf-file' in set(ANNOTATION.keys())
if gtf_index:
    LINK_ANNOTATION    = [os.path.join(PATH_ANNOTATION, 'GTF_Link.gtf')]
    PREPARE_ANNOTATION = [os.path.join(PATH_ANNOTATION, 'Processed_Annotation.rds')]

    ANNOTATE_PEAKS     = expand(os.path.join(PATH_RDS_TEMP,'{name}','{name}.Annotate_Peaks.rds'), name=PEAK_NAMES)
    
#     EXTRACT_SIGNAL_ANNOTATION = expand(os.path.join(PATH_RDS_TEMP,'{name}','{name}.Extract_Signal_Annotation.rds'), name=PEAK_NAMES)    

#     include: os.path.join(RULES_PATH, 'Extract_Signal_Annotation.py')
    include: os.path.join(RULES_PATH, 'Prepare_Annotation.py')
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
# extracts ChIP/Cont signal profiles around the peaks
# EXTRACT_SIGNAL_PEAKS = expand(os.path.join(PATH_RDS_TEMP,'{name}','{name}.Extract_Signal_Peaks.rds'), name=PEAK_NAMES)
# COMMAND = COMMAND + EXTRACT_SIGNAL_PEAKS


# ----------------------------------------------------------------------------- #
# if 'feature_combination' in set(config.keys()) and 'gtf' in set(config['annotation'].keys()


# ----------------------------------------------------------------------------- #
rule all:
    input:
        COMMAND
