# ---------------------------------------------------------------------------- #
import glob
import fnmatch
import os
import re
import sys
import yaml
import xlrd
import csv
import inspect

include: os.path.join(config['locations']['pkglibexecdir'], 'scripts/SnakeFunctions.py')
include: os.path.join(config['locations']['pkglibexecdir'], 'scripts/Check_Config.py')

# ---------------------------------------------------------------------------- #
# check settings and sample_sheet validity
validate_config(config)

localrules: makelinks

SCRIPT_PATH       = os.path.join(config['locations']['pkglibexecdir'], 'scripts/')
RULES_PATH        = os.path.join(config['locations']['pkglibexecdir'], 'Rules/')
SAMPLE_SHEET_FILE = config['locations']['sample-sheet']
REPORT_TEMPLATE   = os.path.join(SCRIPT_PATH,'Sample_Report.rmd')

# ---------------------------------------------------------------------------- #
# read the sample sheet

## Load sample sheet
# either from excel file
if SAMPLE_SHEET_FILE.endswith('.xlsx'):
    with xlrd.open_workbook(file_name) as book:
        # assume that the first book is the sample sheet
        sheet = book.sheet_by_index(0)
        rows = [sheet.row_values(r) for r in range(0, sheet.nrows)]
        header = rows[0]; rows = rows[1:]
        SAMPLE_SHEET = [dict(zip(header, row)) for row in rows]
# or from csv file
elif SAMPLE_SHEET_FILE.endswith('.csv'):
    with open(SAMPLE_SHEET_FILE, 'r') as fp:
        SAMPLE_SHEET = [row for row in csv.DictReader(fp, skipinitialspace=True)]
else:
    raise InputError('File format of the sample_sheet has to be: csv or excel table (xls,xlsx).')


# Convenience function to access fields of sample sheet columns that
# match the predicate.  The predicate may be a string.
def lookup(column, predicate, fields=[]):
    if inspect.isfunction(predicate):
        records = [line for line in SAMPLE_SHEET if predicate(line[column])]
    else:
        records = [line for line in SAMPLE_SHEET if line[column]==predicate]
    return [record[field] for record in records for field in fields]





# ---------------------------------------------------------------------------- #
# Software executables
SOFTWARE = config['tools']


# Per sample software parameters:
# Loops throug the sample sheet and extracts per sample software parameters
# Flattens all samples with custom parameters into one dict - sample names must be unique
custom_param_names = sorted(list(set(config.keys()) -
                     set(['locations', 'general',
                          'execution', 'tools',
                          'hub'])))
CUSTOM_PARAMS = dict()
for param_set in custom_param_names:
    for sample_name in config[param_set].keys():
        sample_set = config[param_set][sample_name]
        if isinstance(sample_set, dict) and 'params' in set(sample_set.keys()):
            CUSTOM_PARAMS[sample_name] = sample_set['params']
        else:
            CUSTOM_PARAMS[sample_name] = None

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
NAMES = [line['SampleName'] for line in SAMPLE_SHEET]

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
FASTQ_FILES = flatten([[sample['Read'], sample['Read2']] for sample in SAMPLE_SHEET])
FASTQC_DICT = {}
for fqfile in FASTQ_FILES:
    if not fqfile == '':
        prefix  = fqfile
        prefix  = re.sub('.fq.*'   , '', prefix)
        prefix  = re.sub('.fastq.*', '', prefix)
        fastqc  = os.path.join(PATH_QC, i, prefix + "_fastqc.zip")
        FASTQC_DICT[prefix]  =  {'fastq'  : os.path.join(PATH_FASTQ, fqfile),
                                 'fastqc' : fastqc}
    # ---------------------------------------------------------------------------- #
# RULE ALL
# Default output files from the pipeline


# ---------------------------------------------------------------------------- #
# width extension parameters for annotation construction
DEFALUT_WIDTH_PARAMS = {
'tss_width':            1000,
'tts_width':            1000,
'tss_wide_width':       10000,
'tts_wide_width':       10000,
'tss_body_upstream':    1000,
'tss_body_downstream':  10000,
'tts_body_upstream':    10000,
'tts_body_downstream':  1000,
'splicing_donor_width': 200,
'splicing_accep_width': 200}

# checks whether the width_params are set, if not
# they are set to DEFALUT_WIDTH_PARAMS
if not 'width_params' in set(PARAMS.keys()):
    PARAMS['width_params'] = DEFALUT_WIDTH_PARAMS
else:
    width_params = PARAMS['width_params']
    if(len(width_params.keys())):
        PARAMS['width_params'] = DEFALUT_WIDTH_PARAMS
    else:
        for i in DEFALUT_WIDTH_PARAMS.keys():
            if not i in set(width_params.keys()):
                width_params[i] = DEFALUT_WIDTH_PARAMS[i]

        PARAMS['width_params'] = width_params


# ---------------------------------------------------------------------------- #
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

COMMAND = GENOME_FASTA + INDEX + BOWTIE2 + CHRLEN + TILLING_WINDOWS + NUCLEOTIDE_FREQ+ BW + LINKS + FASTQC + MULTIQC + BOWTIE2_STATS

# ---------------------------------------------------------------------------- #
# include rules
include: os.path.join(RULES_PATH, 'Mapping.py')
include: os.path.join(RULES_PATH, 'Parse_Bowtie2log.py')
include: os.path.join(RULES_PATH, 'FastQC.py')
include: os.path.join(RULES_PATH, 'MultiQC.py')
include: os.path.join(RULES_PATH, 'BamToBigWig.py')

# ---------------------------------------------------------------------------- #
# Formats the annotation + extracts signal profiles around pre-specified annotation regions
include: os.path.join(RULES_PATH, 'Prepare_Annotation.py')
include: os.path.join(RULES_PATH, 'Extract_Signal_Annotation.py')

LINK_ANNOTATION           = [os.path.join(PATH_ANNOTATION, 'GTF_Link.gtf')]
PREPARE_ANNOTATION        = [os.path.join(PATH_ANNOTATION, 'Processed_Annotation.rds')]
EXTRACT_SIGNAL_ANNOTATION = expand(os.path.join(PATH_RDS_TEMP,'{name}','{name}.Extract_Signal_Annotation.rds'), name=NAMES)

COMMAND = COMMAND + LINK_ANNOTATION + PREPARE_ANNOTATION + EXTRACT_SIGNAL_ANNOTATION

# ---------------------------------------------------------------------------- #
# does the chipqc
include: os.path.join(RULES_PATH, 'ChIPQC.py')
COMMAND = COMMAND + ChIPQC

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# CONDITIONAL OUTPUT FILES

peak_index = 'peak_calling' in set(config.keys())
if peak_index:
    if len(config['peak_calling'].keys()) > 0:
        PEAK_NAMES = config['peak_calling'].keys()
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
if 'idr' in set(config.keys()):
    if len(config['idr'].keys()) > 0:
        IDR = []
        for name in config['idr'].keys():
            IDR = IDR + [os.path.join(PATH_IDR, name, name + ".bed")]
            PEAK_NAME_LIST[name] = IDR[-1]

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


# ---------------------------------------------------------------------------- #
gtf_index = type(ANNOTATION) is str
if peak_index:
    ANNOTATE_PEAKS     = expand(os.path.join(PATH_RDS_TEMP,'{name}','{name}.Annotate_Peaks.rds'), name=PEAK_NAMES)
    include: os.path.join(RULES_PATH, 'Annotate_Peaks.py')
    COMMAND = COMMAND + LINK_ANNOTATION + PREPARE_ANNOTATION + ANNOTATE_PEAKS

# ---------------------------------------------------------------------------- #
if 'feature_combination' in set(config.keys()):
    FEATURE_NAMES = config['feature_combination'].keys()
    if len(FEATURE_NAMES) > 0:

        FEATURE = expand(os.path.join(PATH_RDS_FEATURE,'{name}_FeatureCombination.rds'),
        name = FEATURE_NAMES)

    include: os.path.join(RULES_PATH, 'Feature_Combination.py')
    COMMAND = COMMAND + FEATURE

# ----------------------------------------------------------------------------- #
# REPORT INPUT
SUMMARIZED_DATA_FOR_REPORT = [os.path.join(PATH_ANALYSIS, 'Summarized_Data_For_Report.RDS')]
REPORT_CHUNKS  = {'EXTRACT_SIGNAL_ANNOTATION':'Extract_Signal_Annotation','PEAK_STATISTICS':'Peak_Statistics','ANNOTATE_PEAKS':'Annotate_Peaks','ChIPQC':'ChIPQC'}
REPORT_INPUT   = []
ANALISYS_NAMES = []
for i in REPORT_CHUNKS.keys():
    if i in globals().keys():
        REPORT_INPUT   = REPORT_INPUT   + globals()[i]
        ANALISYS_NAMES = ANALISYS_NAMES + [i]


REPORT         = [os.path.join(PATH_REPORTS, 'ChIP_Seq_Report.html')]
# REPORT_CHUNKS  = ['EXTRACT_SIGNAL_ANNOTATION','PEAK_STATISTICS','ANNOTATE_PEAKS','ChIPQC']
# This lines convert the analysis code chunks from SNAKEMAKE rule language to R language

ANALISYS_NAMES = [REPORT_CHUNKS[i] for i in ANALISYS_NAMES]
include: os.path.join(RULES_PATH, 'Summarize_Data_For_Report.py')
include: os.path.join(RULES_PATH, 'Knit_Report.py')
COMMAND = COMMAND + SUMMARIZED_DATA_FOR_REPORT + REPORT
# ----------------------------------------------------------------------------- #
rule all:
    input:
        COMMAND
