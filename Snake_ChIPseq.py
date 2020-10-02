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

# ---------------------------------------------------------------------------- #
# GLOBAL VARIABLES:
# IMPORTANT defines the compulsory column names for the sample sheet
# if you want to add the compulsory columns for the sample sheet, do it here
STRUCTURE_VARIABLES = {
# Spike-in is a non-compulsory column: values: No/Yes
'SAMPLE_SHEET_COLUMN_NAMES' : ['SampleName', 'Read', 'Read2'],

# defines the allowed execution parameters list for the config file
'SETTING_SUBSECTIONS'       : ['locations', 'general', 'execution', 'tools', 'peak_calling', 'idr', 'hub', 'feature_combination', 'differential_analysis'],

# sets the obligatory files for the pipeline
'OBLIGATORY_FILES'          : ['genome-file','gff-file'],

# obligatory names for report chunks
'REPORT_CHUNKS' : {'EXTRACT_SIGNAL_ANNOTATION':'Extract_Signal_Annotation','PEAK_STATISTICS':'Peak_Statistics','ANNOTATE_PEAKS':'Annotate_Peaks','ChIPQC':'ChIPQC'},

# Types of analysis which are supported
# used in check config
# used to setup global variables (SHIFT)
'ANALYSIS_TYPES' : set(['CHIP','ATAC'])
}


# ---------------------------------------------------------------------------- #
include: os.path.join(config['locations']['pkglibexecdir'], 'scripts/SnakeFunctions.py')
include: os.path.join(config['locations']['pkglibexecdir'], 'scripts/Check_Config.py')
localrules: makelinks
# ---------------------------------------------------------------------------- #
# reads in the sample sheet
# SAMPLE_SHEET is a hardcoded global variable name - it can not change
# check settings and sample_sheet validity
validate_config(config, STRUCTURE_VARIABLES)

SAMPLE_SHEET = read_SAMPLE_SHEET(config)
# ---------------------------------------------------------------------------- #
SCRIPT_PATH       = os.path.join(config['locations']['pkglibexecdir'], 'scripts/')
RULES_PATH        = os.path.join(config['locations']['pkglibexecdir'], 'Rules/')
REPORT_TEMPLATE   = os.path.join(SCRIPT_PATH,'Sample_Report.rmd')
LIB_TYPE          = dict(zip([i['SampleName'] for i in SAMPLE_SHEET],[i['library_type'] for i in SAMPLE_SHEET]))


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

# type of the analysis to run
ANALYSIS     = config['general']['analysis'].upper()

# Sample name definition
PEAK_NAMES   = []
NAMES = [line['SampleName'] for line in SAMPLE_SHEET]

# Directory structure definition
OUTPUT_DIR      = config['locations']['output-dir']
# PATH_FASTQ      = os.path.join(OUTPUT_DIR, 'Fastq')
PATH_TRIMMED    = os.path.join(OUTPUT_DIR, 'Trimmed/Trim_Galore')
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
PATH_RDS_COUNTS     = os.path.join(PATH_RDS, 'Feature_Counts')
PATH_RDS_CHIPQC     = os.path.join(PATH_RDS, 'ChIPQC')
PATH_RDS_TEMP       = os.path.join(PATH_RDS, 'Temp')

# ---------------------------------------------------------------------------- #
# Hardcoded parameters
# determines the shift window for the ChIPQC cross-correlation analysis
# ChIP: 400; ATAC: 2000
if ANALYSIS == 'CHIP':
    SHIFT_WINDOW = 400
    DISCARD_CHRM = PARAMS['chipqc']['discard_chrM']
    ANALYSIS_MODE = "ChIP"

if ANALYSIS == 'ATAC':
    SHIFT_WINDOW = 2000
    DISCARD_CHRM = 'yes'
    ANALYSIS_MODE = "ATAC"
    
    

# ---------------------------------------------------------------------------- #
# ACCESSORY VARIABLES
# Collects the locations of all peaks
PEAK_NAME_LIST = {}

# a list which contains all mapped files
# currently used in Parse_Bowtie2log.py
BAMFILES_LIST = []

# this hash contains the location of the genome link, genome prefix, and bowtie2
# index files for the main and spike-in genomes
# HASH Structure:
# GENOME_HASH[genome_type]['genome_location']
# GENOME_HASH[genome_type]['genome_name']
# GENOME_HASH[genome_type]['genome_prefix']
# GENOME_HASH[genome_type]['genome_link']
# GENOME_HASH[genome_type]['bowtie_index']
GENOME_HASH  = {}

#  contains the genome types for the main and spike-in genome - can be extended
GENOME_TYPES = {'Main' : 'Main','Spike-in' : 'Spike-in'}

# wildcard constraint used for genome building and mapping
GENOME_TYPES_CONSTRAINT = "|".join(list(GENOME_TYPES.values()))

# names for markdown chunks which will be knit
# current names: ChIPQC, Extract_Signal_Annotation, Peak_Statistics, Annotate_Peaks

# Hub variables which describe the types of files that can be used in the hub
TRACK_PATHS = {
    'bigWig' : {'path': os.path.join(PATH_MAPPED, GENOME_TYPES['Main']) , 'suffix': 'bw', 'type':'bigWig'},
    'macs' :{'path': PATH_PEAK, 'suffix': 'bb', 'type':'bigBed'},
    'idr' : {'path': PATH_IDR, 'suffix': 'IDR.narrowPeak'}
}

# ---------------------------------------------------------------------------- #
# Constructs the genome index prefix name
# fills in the GENOME_HASH with paths for the main genome
generate_genome_files(GENOME_ORIG, PATH_INDEX, GENOME_TYPES['Main'], GENOME)
GENOME_MAIN_PREFIX_PATH = GENOME_HASH[GENOME_TYPES['Main']]['genome_prefix']
PATH_MAPPED_MAIN = os.path.join(PATH_MAPPED, GENOME_TYPES['Main'])

# ---------------------------------------------------------------------------- #
# due to the iditotic namig scheme in FASTQC the next lines construct
# FASTQC output files
FASTQC_DICT = {}
for i in NAMES:
    for fqfile in lookup('SampleName',i,['Read','Read2']):
        if not fqfile == '':
            prefix  = fqfile
            prefix  = re.sub('.fq.*'   , '', prefix)
            prefix  = re.sub('.fastq.*', '', prefix)
            fastqc  = os.path.join(PATH_QC, i, prefix + "_fastqc.zip")
            FASTQC_DICT[prefix]  =  {'fastq'  : os.path.join(PATH_FASTQ, fqfile),
                                     'fastqc' : fastqc}

# ---------------------------------------------------------------------------- #
# due to the different names for trimmmed output files next lines construct
# Trim Galore output files
TRIM_GALORE_DICT = {}
TRIM_GALORE_FILES = {}
for name in NAMES:
    TRIM_GALORE_DICT[name] = get_trimming_dict(name)
    TRIM_GALORE_FILES[name] = flatten([TRIM_GALORE_DICT[name][rep]['trimmed'] for rep in TRIM_GALORE_DICT[name].keys()])




# ---------------------------------------------------------------------------- #
# width extension parameters for annotation construction
DEFAULT_WIDTH_PARAMS = {
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
# they are set to DEFAULT_WIDTH_PARAMS
if not 'width_params' in set(PARAMS.keys()):
    PARAMS['width_params'] = DEFAULT_WIDTH_PARAMS
else:
    width_params = PARAMS['width_params']
    if(len(width_params.keys())):
        PARAMS['width_params'] = DEFAULT_WIDTH_PARAMS
    else:
        for i in DEFAULT_WIDTH_PARAMS.keys():
            if not i in set(width_params.keys()):
                width_params[i] = DEFAULT_WIDTH_PARAMS[i]

        PARAMS['width_params'] = width_params


# ---------------------------------------------------------------------------- #
# check wether quality filtering should be performed, if yes
# change the suffix of the output bam file

BAM_SUFFIX = ''

if PARAMS['bam_filter']['mapq']:
    if PARAMS['bam_filter']['mapq'] > 0:
        BAM_SUFFIX = '.q{}'.format(PARAMS['bam_filter']['mapq'])
        #BAM_SUFFIX = '.qfilter'

# check wether deduplication should be performed, if yes
# change the suffix of the output bam file
if PARAMS['bam_filter']['deduplicate']:
    BAM_SUFFIX = BAM_SUFFIX + ".deduplicated.sorted.bam"
else:
    BAM_SUFFIX = BAM_SUFFIX + ".sorted.bam"


# ---------------------------------------------------------------------------- #
# Inline definition and description of targets
targets = {
    # rule to print all rule descriptions
    'help': {
        'description': "Print all rules and their descriptions.",
        'files': []
    }
}

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# GENERAL MAPPING OUTPUT FILES
GENOME_FASTA    = [GENOME_HASH[GENOME_TYPES['Main']]['genome_link']]
INDEX           = [GENOME_HASH[GENOME_TYPES['Main']]['bowtie_index']]
TRIMMING        = flatten(TRIM_GALORE_FILES.values())

targets['trimming'] = {
        'description': "Trim reads with Trim Galore.",
        'files':
            TRIMMING
}

# maps the reads to the main genome
BOWTIE2         = expand(
    os.path.join(PATH_MAPPED, GENOME_TYPES['Main'],
    "{name}", "{name}" + BAM_SUFFIX + ".bai"), name=NAMES)
    
BAMFILES_LIST.append(BOWTIE2)

BOWTIE2_STATS   = [os.path.join(PATH_RDS, "BowtieLog.rds")]
CHRLEN          = [GENOME_MAIN_PREFIX_PATH + '.chrlen.txt']
TILLING_WINDOWS = [GENOME_MAIN_PREFIX_PATH + '.GenomicWindows.GRanges.rds']
NUCLEOTIDE_FREQ = [GENOME_MAIN_PREFIX_PATH + '.NucleotideFrequency.GRanges.rds']
FASTQC          = [FASTQC_DICT[i]['fastqc'] for i in list(FASTQC_DICT.keys())]
MULTIQC         = [os.path.join(PATH_REPORTS, "multiqc.html")]
ChIPQC          = expand(
    os.path.join(
        PATH_RDS_CHIPQC, 
        GENOME_TYPES['Main'], 
        GENOME_HASH['Main']['genome_name'], 
        "{name}_ChIPQC.rds"), 
    name=NAMES)
BW              = expand(
    os.path.join(PATH_MAPPED, GENOME_TYPES['Main'],
    "{name}", "{name}.bw"),  name=NAMES)
LINKS           = expand(
    os.path.join(PATH_BW, GENOME_TYPES['Main'], "{ex_name}.bw"),  ex_name=NAMES)



targets['mapping'] = {
        'description': "Produce the bowtie2 mapping results in BAM format.",
        'files':
            BOWTIE2
}

targets['export-bw'] = {
        'description': "Take the bowtie2 mapping results in BAM format and create bigWig Tracks.",
        'files':
            BW + LINKS
}

targets['multiqc'] = {
        'description': "Get multiQC report based on bowtie2 alignments and fastQC reports.",
        'files':
            FASTQC +  MULTIQC
}

# defines the rules for spike in quantification
if 'spikein-file' in set(config['locations'].keys()):

    # defines the spiked samples
    NAMES_SPIKEIN = []
    for line in SAMPLE_SHEET:
        if 'Spike-in' in set(line.keys()):
            if line['Spike-in'].upper() == 'YES':
                NAMES_SPIKEIN.append(line['SampleName'])

    if(len(NAMES_SPIKEIN) > 0):

        GENOME_SPIKEIN  = config['locations']['spikein-file']

        # generates the spike in name, if the name is not supplied
        if 'spikein_name' in set(config['general'].keys()):
            GENOME_NAME_SPIKEIN  = config['general']['spikein_name']
        else:
            GENOME_NAME_SPIKEIN  = os.path.basename(GENOME_SPIKEIN)

        # fills in the genome hash with the spike-in genome
        generate_genome_files(GENOME_SPIKEIN, PATH_INDEX, GENOME_TYPES['Spike-in'], GENOME_NAME_SPIKEIN)

        # maps the reads to the spike in genome
        BOWTIE2_SPIKEIN  = expand(os.path.join(PATH_MAPPED, GENOME_TYPES['Spike-in'], "{name}", "{name}" + BAM_SUFFIX + ".bai"), name=NAMES_SPIKEIN)

        BAMFILES_LIST.append(BOWTIE2_SPIKEIN)
        SPIKE_IN_RULES = [
            GENOME_HASH[GENOME_TYPES['Spike-in']]['genome_link'],
            GENOME_HASH[GENOME_TYPES['Spike-in']]['bowtie_index'],
        ] + BOWTIE2_SPIKEIN

        targets['spike-in'] = {
            'description' : 'map reads to the spike-in genome and calculate normalization factors',
            'files':
                SPIKE_IN_RULES + BW + LINKS
        }



# ---------------------------------------------------------------------------- #
# include rules
include: os.path.join(RULES_PATH, 'Trimming.py')
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

# ---------------------------------------------------------------------------- #
# does the chipqc
include: os.path.join(RULES_PATH, 'ChIPQC.py')

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


        targets['peak-calling'] = {
            'description': "Perform peak calling based on peak_calling section.",
            'files':
             MACS +  QSORT +  PEAK_STATISTICS
            }


# # ----------------------------------------------------------------------------- #
if 'idr' in set(config.keys()):
    if len(config['idr'].keys()) > 0:
        IDR = []
        for name in config['idr'].keys():
            IDR = IDR + [os.path.join(PATH_IDR, name, name + ".bed")]
            PEAK_NAME_LIST[name] = IDR[-1]

        include: os.path.join(RULES_PATH, 'IDR.py')

        targets['idr'] = {
            'description': "Control reproducibilty of peak calling based on idr section.",
            'files':IDR
            }

# # ----------------------------------------------------------------------------- #
HUB_NAME = None
if 'hub' in set(config.keys()):
    HUB_NAME = config['hub']['name']
    HUB = [os.path.join(PATH_HUB, HUB_NAME, 'done.txt')]
    BB  = expand(os.path.join(PATH_PEAK,  "{name}", "{name}.bb"),  name=config['peak_calling'].keys())

    include: os.path.join(RULES_PATH, 'UCSC_Hub.py')

    targets['hub'] = {
            'description': "Generate UCSC track hub based on tracks defined in hub section.",
            'files': BB + HUB
            }


# ---------------------------------------------------------------------------- #
gtf_index = type(ANNOTATION) is str
if peak_index:
    ANNOTATE_PEAKS     = expand(os.path.join(PATH_RDS_TEMP,'{name}','{name}.Annotate_Peaks.rds'), name=PEAK_NAMES)
    include: os.path.join(RULES_PATH, 'Annotate_Peaks.py')

# ---------------------------------------------------------------------------- #
if 'feature_combination' in set(config.keys()):
    FEATURE_NAMES = config['feature_combination'].keys()
    if len(FEATURE_NAMES) > 0:
        FEATURE = expand(os.path.join(PATH_RDS_FEATURE,'{name}_FeatureCombination.rds'),
                         name = FEATURE_NAMES)

        include: os.path.join(RULES_PATH, 'Feature_Combination.py')

        targets['feature-combination'] = {
            'description': "Identify overlapping features based on feature_combination section.",
            'files': FEATURE
            }

# ---------------------------------------------------------------------------- #
if 'differential_analysis' in set(config.keys()):
    DIFF_ANALYSIS_NAMES = config['differential_analysis'].keys()
    if len(DIFF_ANALYSIS_NAMES) > 0:

        include: os.path.join(RULES_PATH, 'Differential_Analysis.py')

        DIFF_ANALYSIS_CONSENSENSUS = []
        DIFF_ANALYSIS_COUNTS = []

        for diffAnn in DIFF_ANALYSIS_NAMES:

            diffAnnDict = config['differential_analysis'][diffAnn]

            # if no peaks given call joint peaks for given samples  
            if ( not 'Peakset' in diffAnnDict ) or ( not diffAnnDict['Peakset']):
                DIFF_ANALYSIS_CONSENSENSUS += [ 
                        os.path.join(PATH_PEAK,  diffAnn, diffAnn + "_peaks.narrowPeak"),
                        os.path.join(PATH_PEAK,  diffAnn, diffAnn + "_qsort.bed" )
                ]
                ## TODO find less hacky way to add to peak dict 
                config['peak_calling'].update(
                        { diffAnn : {
                    'ChIP': diffAnnDict['Case'] + diffAnnDict['Control'],
                    'Cont': None}})
            else:
                DIFF_ANALYSIS_CONSENSENSUS += expand(os.path.join(PATH_RDS_FEATURE,'{name}_FeatureCombination.{type}'),
                     name = diffAnn, type = ['rds','txt'])
                
                if not 'feature_combination' in set(config.keys()):
                    include: os.path.join(RULES_PATH, 'Feature_Combination.py')
                    config.update({'feature_combination':dict()}) 
                    ## TODO find less hacky way to add to feature_combination dict 
                config['feature_combination'].update(
                        { diffAnn : diffAnnDict['Peakset']})
                print(config['feature_combination'][diffAnn])

                include: os.path.join(RULES_PATH, 'Feature_Combination.py')

        print(DIFF_ANALYSIS_CONSENSENSUS)
        targets['differential-analysis'] = {
            'description': "Identify peaks with differential read occupancy.",
            'files': DIFF_ANALYSIS_CONSENSENSUS
            }


# ----------------------------------------------------------------------------- #
# REPORT INPUT
SUMMARIZED_DATA_FOR_REPORT = [os.path.join(PATH_ANALYSIS, 'Summarized_Data_For_Report.RDS')]
REPORT_CHUNKS  = STRUCTURE_VARIABLES['REPORT_CHUNKS']
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

targets['final-report'] = {
'description': "Produce a comprehensive report.",
'files': SUMMARIZED_DATA_FOR_REPORT + REPORT
}

# ----------------------------------------------------------------------------- #
# COMPLETE EXECUTION
ALL_FILES = list(chain.from_iterable([targets[name]['files'] for name in list(targets.keys())]))
targets['complete'] = {
'description': "Run the full pipeline with all available targetds.  This is the default target.",
'files': ALL_FILES
}
# ----------------------------------------------------------------------------- #
# TARGETTED EXECUTION
# Selected output files from the above set.
selected_targets = config['execution']['target'] or ['complete']

# FIXME: the list of files must be flattened twice(!).  We should make
# sure that the targets really just return simple lists.
from itertools import chain
wrong_target = []
for selection in selected_targets:
    if not selection in targets.keys():
            wrong_target.append(selection)
    if wrong_target:
        sys.exit(''.join(
            ['This is not a supported targed:\n {}\n'.format(wrong_target),
             'Consider one of these:\n {}\n'.format(list(targets.keys()))] )
            )
OUTPUT_FILES = list(chain.from_iterable([targets[name]['files'] for name in selected_targets]))

rule all:
  input: OUTPUT_FILES

rule help:
  run:
    for key in targets.keys():
      print('{}:\n  {}'.format(key, targets[key]['description']))

# Record any existing output files, so that we can detect if they have
# changed.
expected_files = {}
onstart:
    if OUTPUT_FILES:
        for name in OUTPUT_FILES:
            if os.path.exists(name):
                expected_files[name] = os.path.getmtime(name)

# Print generated target files.
onsuccess:
    if OUTPUT_FILES:
        # check if any existing files have been modified
        generated = []
        for name in OUTPUT_FILES:
            if name not in expected_files or os.path.getmtime(name) != expected_files[name]:
                generated.append(name)
        if generated:
            print("The following files have been generated:")
            for name in generated:
                print("  - {}".format(name))
