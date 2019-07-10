
"""
Snakefile for pigx-scrnaseq pipeline
"""

# ----------------------------------------------------------------------------- #
# libraries and constants
import glob
import os
import re
import subprocess
import yaml
import csv
import inspect
import magic as mg
import sys
import pandas as pd

PATH_SCRIPT = os.path.join(config['locations']['pkglibexecdir'], 'scripts')

# loads function for running Rscripts
include: os.path.join(PATH_SCRIPT, 'Run_Rscript.py')
# functions for java input
include: os.path.join(PATH_SCRIPT, 'Accessory_Functions.py')
include: os.path.join(PATH_SCRIPT, 'validate_input.py')
validate_config(config)
# class for SAMPLE_SHEET
include: os.path.join(PATH_SCRIPT, 'Sample_Sheet_Class.py')

# ----------------------------------------------------------------------------- #
# Software parameters
SOFTWARE = config['tools']
# variables

GENOME_NAME_PRIMARY = config['annotation']['primary']['genome']['name']
REFERENCE_NAMES     = [GENOME_NAME_PRIMARY]
COVARIATES          = config['covariates']

# ----------------------------------------------------------------------------- #
# adapter locations
ADAPTER_PARAMETERS = config['adapter_parameters']


# ----------------------------------------------------------------------------- #
# PATHS
OUTPUT_DIR           = config['locations']['output-dir']
PATH_FASTQ           = config['locations']['reads-dir']
TEMPDIR              = config['locations']['tempdir']
PATH_ANNOTATION      = os.path.join(OUTPUT_DIR, 'Annotation')
PATH_FASTQC          = os.path.join(OUTPUT_DIR, 'FASTQC')
PATH_MAPPED          = os.path.join(OUTPUT_DIR, 'Mapped')
PATH_LOG             = os.path.join(OUTPUT_DIR, 'Log')
PATH_SAMPLE_SHEET    = config['locations']['sample-sheet']
PATH_RSCRIPT         = SOFTWARE['Rscript']['executable']

PATH_ANNOTATION_PRIMARY = os.path.join(PATH_ANNOTATION, GENOME_NAME_PRIMARY)
PATH_REFERENCE_PRIMARY  = config['annotation']['primary']['genome']['fasta']
PATH_GTF_PRIMARY        = config['annotation']['primary']['gtf']

GENOME_NAME_MIX = None
if config['annotation']['secondary']:
    GENOME_NAME_SECONDARY    = config['annotation']['secondary']['genome']['name']
    PATH_REFERENCE_SECONDARY = config['annotation']['secondary']['genome']['fasta']
    PATH_GTF_SECONDARY       = config['annotation']['secondary']['gtf']

    GENOME_NAME_MIX = GENOME_NAME_PRIMARY + '_' + GENOME_NAME_SECONDARY
    PATH_ANNOTATION_MIX = os.path.join(PATH_ANNOTATION, GENOME_NAME_MIX)
    REFERENCE_NAMES = REFERENCE_NAMES + [GENOME_NAME_MIX]

## Load sample sheet
SAMPLE_SHEET = experiment(config = config)
SAMPLE_SHEET.init_SAMPLE_SHEET(PATH_SAMPLE_SHEET)

SAMPLE_NAMES = SAMPLE_SHEET.fetch_sample_names()

# ----------------------------------------------------------------------------- #
# sets the temporrary directory to default in the working directory if the tempdir does not exist
if TEMPDIR == None:
    if "TMPDIR" in os.environ.keys():
        TEMPDIR = os.environ['TMPDIR']
    else:
        TEMPDIR = '/tmp'

# ----------------------------------------------------------------------------- #
# RULES

# ----------------------------------------------------------------------------- #
# Link primary reference
# TODO : make extension dynamic (to accept both fasta and fasta.gz files)
LINK_REFERENCE_PRIMARY   = os.path.join(PATH_ANNOTATION_PRIMARY,  GENOME_NAME_PRIMARY + '.fasta')
LINK_GTF_PRIMARY         = os.path.join(PATH_ANNOTATION_PRIMARY, GENOME_NAME_PRIMARY + '.gtf')


# ----------------------------------------------------------------------------- #
# Combine primary and secondary reference genomes
COMBINE_REFERENCE = []
GENOME_SECONDARY_IND = not GENOME_NAME_MIX == None
if GENOME_SECONDARY_IND:
    PATH_REFERENCE_MIX = os.path.join(PATH_ANNOTATION_MIX, GENOME_NAME_MIX + '.fasta')
    PATH_GTF_MIX = os.path.join(PATH_ANNOTATION_MIX, GENOME_NAME_MIX + '.gtf')
    COMBINE_REFERENCE = COMBINE_REFERENCE + [PATH_REFERENCE_MIX, PATH_GTF_MIX]

# ----------------------------------------------------------------------------- # FastQC
FASTQC_prep  = list(SAMPLE_SHEET.SAMPLE_SHEET['reads']) + list(SAMPLE_SHEET.SAMPLE_SHEET['barcode'])
FASTQC = [os.path.join(PATH_FASTQC, file + ".fastqc.done") for file in FASTQC_prep]

# ----------------------------------------------------------------------------- #
# Change reference gene_name to gene_id
PATH_GTF_PRIMARY_ID = expand(os.path.join(PATH_ANNOTATION_PRIMARY, "{name}" + "gene_id.gtf"), name=REFERENCE_NAMES)

# ----------------------------------------------------------------------------- #
# STAR INDEX
MAKE_STAR_INDEX = expand(os.path.join(PATH_ANNOTATION, '{genome}','STAR_INDEX','done.txt'), genome = REFERENCE_NAMES)


# -----------------------------------------------------------------------------
#
# merge reads
# use for loop - expand gives wrong results
MERGE_TECHNICAL_REPLICATES = []
for sample_name in SAMPLE_NAMES:
    print(sample_name)
    MERGE_TECHNICAL_REPLICATES.append(os.path.join(PATH_MAPPED, sample_name, SAMPLE_SHEET.fetch_field(sample_name,'reads_merged')))
    MERGE_TECHNICAL_REPLICATES.append(os.path.join(PATH_MAPPED, sample_name, SAMPLE_SHEET.fetch_field(sample_name,'barcode_merged')))

# ----------------------------------------------------------------------------- #
# MAPPING
MAP_scRNA = expand(os.path.join(PATH_MAPPED, "{name}", "{genome}", "{name}.bam"), genome = REFERENCE_NAMES, name = SAMPLE_NAMES)


# ----------------------------------------------------------------------------- #
# Number of reads per cell calculation
BAM_HISTOGRAM = expand(os.path.join(PATH_MAPPED, "{name}", "{genome}",'{name}_{genome}_BAMTagHistogram.txt'), genome = REFERENCE_NAMES, name = SAMPLE_NAMES)

# ----------------------------------------------------------------------------- #
# Number of reads per cell calculation
FIND_READ_CUTOFF = expand(os.path.join(PATH_MAPPED, "{name}", "{genome}",'{name}_{genome}_ReadCutoff.yaml'), genome = REFERENCE_NAMES, name = SAMPLE_NAMES)


# ----------------------------------------------------------------------------- #
# Reads matrix
READS_MATRIX = expand(os.path.join(PATH_MAPPED, "{name}", "{genome}",'{name}_{genome}_READS.Matrix.txt'), genome = REFERENCE_NAMES, name = SAMPLE_NAMES)

# ----------------------------------------------------------------------------- #
# UMI matrix
UMI = expand(os.path.join(PATH_MAPPED, "{name}", "{genome}",'{name}_{genome}_UMI.Matrix.txt'), genome = REFERENCE_NAMES, name = SAMPLE_NAMES)

# UMI matrix in loom format
UMI_LOOM =  expand(os.path.join(PATH_MAPPED, "{name}", "{genome}",'{name}_{genome}_UMI.Matrix.loom'), genome = REFERENCE_NAMES, name = SAMPLE_NAMES)

# Combined UMI matrices in loom format
COMBINED_UMI_MATRICES = expand(os.path.join(PATH_MAPPED, "{genome}_UMI.loom"), genome = REFERENCE_NAMES)

# ----------------------------------------------------------------------------- #
# READ statistics
READ_STATISTICS = expand(os.path.join(PATH_MAPPED, "{name}", "{genome}",'{name}_{genome}_ReadStatistics.txt'), genome = REFERENCE_NAMES, name = SAMPLE_NAMES)


# ----------------------------------------------------------------------------- #
# DOWNSTREAM statistics
# IMPORTANT - STILL NOT IMPLEMENTED
DOWNSTREAM_STATISTICS = expand(os.path.join(PATH_MAPPED, "{name}", "{genome}",'{name}_{genome}_DownstreamStatistics.txt'), genome = REFERENCE_NAMES, name = SAMPLE_NAMES)


# ----------------------------------------------------------------------------- #
# Import and preprocess the combined loom files and save as SingleCellExperiment.RDS objects.
SCE_RDS_FILES = expand(os.path.join(PATH_MAPPED, "{genome}.SingleCellExperiment.RDS"), genome = REFERENCE_NAMES)


# ----------------------------------------------------------------------------- #
# Seurat RDS files
SEURAT_RDS_FILES = expand(os.path.join(PATH_MAPPED, "{genome}.Seurat.RDS"), genome = REFERENCE_NAMES)

# ----------------------------------------------------------------------------- #
## Using the preprocessed SingleCellExperiment.RDS file, generates a self-contained HTML report
REPORT_FILES = expand(os.path.join(PATH_MAPPED, "{genome}.scRNA-Seq.report.html"), genome = REFERENCE_NAMES)
# ----------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------- #
# Bam To BigWig
BIGWIG = expand(os.path.join(PATH_MAPPED, "{name}", "{genome}",'{name}_{genome}.bw'), genome = REFERENCE_NAMES, name = SAMPLE_NAMES)

# ----------------------------------------------------------------------------- #
RULE_ALL = []
RULE_ALL = RULE_ALL + [LINK_REFERENCE_PRIMARY, LINK_GTF_PRIMARY]

if len(COMBINE_REFERENCE) > 0:
    RULE_ALL = RULE_ALL + COMBINE_REFERENCE

#RULE_ALL = RULE_ALL + DICT + REFFLAT + MAKE_STAR_INDEX + FASTQC + MERGE_FASTQ_TO_BAM + MERGE_BAM_PER_SAMPLE + MAP_scRNA + BAM_HISTOGRAM + FIND_READ_CUTOFF + READS_MATRIX + UMI + READ_STATISTICS  + BIGWIG + UMI_LOOM + COMBINED_UMI_MATRICES + SCE_RDS_FILES + SEURAT_RDS_FILES + REPORT_FILES

RULE_ALL = RULE_ALL + MAKE_STAR_INDEX + MERGE_TECHNICAL_REPLICATES


# ----------------------------------------------------------------------------- #
rule all:
    input:
        RULE_ALL

# ----------------------------------------------------------------------------- #
# links the primary annotation to the ./Annotation folder
rule link_primary_annotation:
    input:
        gtf   = PATH_GTF_PRIMARY,
        fasta = PATH_REFERENCE_PRIMARY
    output:
        gtf   = LINK_GTF_PRIMARY,
        fasta = LINK_REFERENCE_PRIMARY
    params:
        threads = config['execution']['rules']['link_primary_annotation']['threads'],
        mem     = config['execution']['rules']['link_primary_annotation']['memory'],
        #ln      = SOFTWARE['ln']['executable']
        ln = 'ln'
    message:
        """
            Linking primary reference files:
                gtf:
                    file: {input.gtf}
                    link: {output.gtf}
                fasta:
                    file: {input.fasta}
                    link: {output.fasta}
        """
    shell:"""
        {params.ln} -s {input.gtf} {output.gtf}
        {params.ln} -s {input.fasta} {output.fasta}
    """



# ----------------------------------------------------------------------------- #
if GENOME_SECONDARY_IND:
    rule combine_reference:
        input:
            primary   =  LINK_REFERENCE_PRIMARY,
            secondary =  PATH_REFERENCE_SECONDARY
        output:
            outfile = PATH_REFERENCE_MIX
        params:
            threads = config['execution']['rules']['combine_reference']['threads'],
            mem     = config['execution']['rules']['combine_reference']['memory'],
            genome_name_primary   = GENOME_NAME_PRIMARY,
            genome_name_secondary = GENOME_NAME_SECONDARY,
            perl = SOFTWARE['perl']['executable'],
            perl_args = SOFTWARE['perl']['executable']['args']
        message:
            """
                Combining fasta files:
                    primary   : {input.primary}
                    secondary : {input.secondary}
                    output: {output.outfile}
            """
        shell:"""
            {params.perl} {params.perl_args} -pe 's|^>|>{params.genome_name_primary}|' < {input.primary} > {output.outfile}
            {params.perl} {params.perl_args} -pe 's|^>|>{params.genome_name_secondary}|' < {input.secondary} >> {output.outfile}
    """

# ----------------------------------------------------------------------------- #
# changes the gene_name field in the GTF file to the gene_id
# this is required for droptools counting
rule change_gtf_id:
    input:
        infile = os.path.join(PATH_ANNOTATION, '{genome}', '{genome}.gtf')
    output:
        outfile = os.path.join(PATH_ANNOTATION, '{genome}', '{genome}.gene_id.gtf')
    params:
        threads = config['execution']['rules']['change_gtf_id']['threads'],
        mem     = config['execution']['rules']['change_gtf_id']['memory'],
        script  = PATH_SCRIPT,
        Rscript = PATH_RSCRIPT,
        perl    = SOFTWARE['perl']['executable']
    message:
        """
            Changing GTF id:
                input  : {input}
                output : {output}
        """
    run:
        RunRscript(input, output, params, params.script, 'change_gtf_id.R')

        # removes header from the gtf file
        command = " ".join([
            params.perl, "-pi -e 's|^#.+$||'", output.outfile
        ])
        shell(command)


# ----------------------------------------------------------------------------- #
# STAR INDEX
###
rule make_star_reference:
    input:
        fasta = os.path.join(PATH_ANNOTATION, '{genome}', '{genome}.fasta'),
        gtf   = rules.change_gtf_id.output.outfile
        # gtf   = os.path.join(PATH_ANNOTATION, '{genome}', '{genome}.gtf')
    output:
        outfile = os.path.join(PATH_ANNOTATION, '{genome}','STAR_INDEX','done.txt')
    params:
        outdir  = os.path.join(PATH_ANNOTATION, '{genome}','STAR_INDEX'),
        star    = SOFTWARE['star']['executable'],
        star_args = SOFTWARE['star']['args'],
        threads = config['execution']['rules']['make_star_reference']['threads'],
        mem     = config['execution']['rules']['make_star_reference']['memory']

    log:
        logfile = os.path.join(PATH_LOG, '{genome}.make_star_reference.log')
    message:"""
        Star reference:
            input:
                fasta : {input.fasta}
                gtf   : {input.gtf}
        """
    run:
        # star genome build command
        command = " ".join(
        [params.star, params.star_args,
         '--runMode genomeGenerate',
        '--genomeDir',        str(params.outdir),
        '--runThreadN',       str(params.threads),
        '--genomeFastaFiles', str(input.fasta),
        '--sjdbGTFfile',      str(input.gtf),
        '--sjdbOverhang',     '99',
        '2>',log.logfile
        ])

        # touch command for the star index
        command_touch = " ".join([
        'touch', str(output.outfile),
        '2>',    str(log.logfile)
        ])

        command_final = command + ';' + command_touch
        print(command_final, file=sys.stderr)
        shell(command_final)

# ----------------------------------------------------------------------------- #
# GIVEN PRIMARY AND SECONDARY GTF, COMBINES THEM INTO ONE GTF FILE
if GENOME_SECONDARY_IND:
    rule combine_gtf:
        input:
            primary   =  LINK_GTF_PRIMARY,
            secondary =  PATH_GTF_SECONDARY
        output:
            outfile = PATH_GTF_MIX
        params:
            threads = config['execution']['rules']['combine_gtf']['threads'],
            mem     = config['execution']['rules']['combine_gtf']['memory'],
            genome_name_primary   = GENOME_NAME_PRIMARY,
            genome_name_secondary = GENOME_NAME_SECONDARY,
            perl = SOFTWARE['perl']['executable'],
            perl_args = SOFTWARE['perl']['executable']['args']
        message:
            """
                Combining gtf files:
                    primary   : {input.primary}
                    secondary : {input.secondary}
                    output: {output.outfile}
            """
        shell:"""
            {params.perl} {params.perl_args} -pe 's|^|{params.genome_name_primary}|' < {input.primary} > {output.outfile}
            {params.perl} {params.perl_args} -pe 's|^|{params.genome_name_secondary}|' < {input.secondary} >> {output.outfile}
    """

# ----------------------------------------------------------------------------- #
def fetch_fastq_technical_replicates(wc):
    sample_name = str(wc.name)
    if(str(wc.type) == 'R1'):
        input_fastq = SAMPLE_SHEET.fetch_barcode_path(sample_name)

    if(str(wc.type) == 'R2'):
        input_fastq = SAMPLE_SHEET.fetch_reads_path(sample_name)

    if not input_fastq.__class__ == list:
        input_fastq = list(input_fastq)
    return(input_fastq)

rule merge_bam_per_sample:
    input:
        bam_files = fetch_all_bams
    output:
        outfile = os.path.join(PATH_MAPPED, "{name}","{name}" + ".Merged.bam")
    log:
        log = os.path.join(PATH_LOG, '{name}.merge_bam_per_sample.log')
    params:
        samtools      = SOFTWARE['samtools']['executable'],
        samtools_args = SOFTWARE['samtools']['args']
    run:
        in_files = ' '.join(input.bam_files)
        command = ' '.join([
            params.samtools,
            params.samtools_args,
            'cat -o',
            output.outfile,
            in_files,
            '2>' + str(log.log)
            ])
        print(command, file=sys.stderr)
        shell(command)

rule merge_technical_replicates:
    input:
        infiles = fetch_fastq_technical_replicates
    output:
        outfile = os.path.join(PATH_MAPPED, "{name}", ("{name}" + "_{type,(R1)|(R2)}" + ".fastq.gz"))
    params:
        threads = config['execution']['rules']['merge_technical_replicates']['threads'],
        mem     = config['execution']['rules']['merge_technical_replicates']['memory'],
        tempdir = TEMPDIR,
        cat     = SOFTWARE['cat']['executable']
    log:
        log = os.path.join(PATH_LOG, '{name}.{type,(R1)|(R2)}.merge_technical_replicates.log')
    message:"""
            merge_technical_replicates:
                input:   {input}
                output : {output}
        """
    run:

        command = ' '.join([
            params.cat,
            " ".join(input.infiles),
            '>'  + str(output.outfile),
            '2>' + str(log.log)
        ])
        print(command, file=sys.stderr)
        shell(command)
