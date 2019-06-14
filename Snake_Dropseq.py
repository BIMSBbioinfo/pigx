
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
REFERENCE_NAMES = [GENOME_NAME_PRIMARY]
COVARIATES = config['covariates']

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

SAMPLE_NAMES = list(set(SAMPLE_SHEET.list_attr('sample_name')))

# ----------------------------------------------------------------------------- #
# check for compatible technology methods
methods = set(ADAPTER_PARAMETERS.keys())
for method in set(SAMPLE_SHEET.list_attr('method')):
    if not method in methods:
        message = 'Sample sheet contains unknown method:' + method + '\n'
        message = message + 'Supported methods are:' + " ".join(list(methods)) + '\n'
        sys.exit(message)

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


# ----------------------------------------------------------------------------- #
# REFFLAT and DICT
REFFLAT = [os.path.join(PATH_ANNOTATION_PRIMARY, GENOME_NAME_PRIMARY + '.refFlat')]
DICT    = [os.path.join(PATH_ANNOTATION_PRIMARY, GENOME_NAME_PRIMARY + '.dict')]

if GENOME_SECONDARY_IND:
    REFFLAT = REFFLAT + [os.path.join(PATH_ANNOTATION_MIX, GENOME_NAME_MIX + '.refFlat')]
    DICT    = DICT    + [os.path.join(PATH_ANNOTATION_MIX, GENOME_NAME_MIX + '.dict')]


# ----------------------------------------------------------------------------- # FastQC
FASTQC_prep  = SAMPLE_SHEET.list_attr('reads') + SAMPLE_SHEET.list_attr('barcode')
FASTQC = [os.path.join(PATH_FASTQC, file + ".fastqc.done") for file in FASTQC_prep]

# ----------------------------------------------------------------------------- #
# Change reference gene_name to gene_id
PATH_GTF_PRIMARY_ID = expand(os.path.join(PATH_ANNOTATION_PRIMARY, "{name}" + "gene_id.gtf"), name=REFERENCE_NAMES)

# ----------------------------------------------------------------------------- #
# STAR INDEX
MAKE_STAR_INDEX = expand(os.path.join(PATH_ANNOTATION, '{genome}','STAR_INDEX','done.txt'), genome = REFERENCE_NAMES)

# ----------------------------------------------------------------------------- #
# MERGE BARODE AND READS FASTQ FILES
SAMPLE_SHEET.merged_bam(PATH_MAPPED)
MERGE_FASTQ_TO_BAM = SAMPLE_SHEET.list_attr('Merged_BAM')

# ----------------------------------------------------------------------------- #
# MERGE ALL BAM FILES PER SAMPLE
MERGE_BAM_PER_SAMPLE = expand(os.path.join(PATH_MAPPED, "{name}", "{name}" + ".Merged.bam"), name = SAMPLE_NAMES)

# ----------------------------------------------------------------------------- #
# MAPPING
MAP_scRNA = expand(os.path.join(PATH_MAPPED, "{name}", "{genome}", "star_gene_exon_tagged.bam"), genome = REFERENCE_NAMES, name = SAMPLE_NAMES)


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

RULE_ALL = RULE_ALL + DICT + REFFLAT + MAKE_STAR_INDEX + FASTQC + MERGE_FASTQ_TO_BAM + MERGE_BAM_PER_SAMPLE + MAP_scRNA + BAM_HISTOGRAM + FIND_READ_CUTOFF + READS_MATRIX + UMI + READ_STATISTICS  + BIGWIG + UMI_LOOM + COMBINED_UMI_MATRICES + SCE_RDS_FILES + REPORT_FILES

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
        mem     = config['execution']['rules']['link_primary_annotation']['memory']
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
        ln -s {input.gtf} {output.gtf}
        ln -s {input.fasta} {output.fasta}
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
# STAR INDEX
rule make_star_reference:
    input:
        fasta = os.path.join(PATH_ANNOTATION, '{genome}', '{genome}.fasta'),
        gtf   = os.path.join(PATH_ANNOTATION, '{genome}', '{genome}.gtf'),
    output:
        outfile = os.path.join(PATH_ANNOTATION, '{genome}','STAR_INDEX','done.txt')
    params:
        outdir  = os.path.join(PATH_ANNOTATION, '{genome}','STAR_INDEX'),
        star    = SOFTWARE['star']['executable'],
        star_args = SOFTWARE['star']['args'],
        threads = config['execution']['rules']['make_star_reference']['threads'],
        mem     = config['execution']['rules']['make_star_reference']['memory']
    log:
        os.path.join(PATH_LOG, '{genome}.make_star_reference.log')
    message:"""
        Star reference:
            input:
                fasta : {input.fasta}
                gtf   : {input.gtf}
        """
    shell:"""
        {params.star} {params.star_args} --runMode genomeGenerate --genomeDir {params.outdir} --genomeFastaFiles {input.fasta} --runThreadN {params.threads} --sjdbGTFfile {input.gtf} --sjdbOverhang 99
        touch {output.outfile} 2> {log}
"""

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

rule fasta_dict:
    input:
        infile = os.path.join(PATH_ANNOTATION, '{genome}', '{genome}.fasta')
    output:
        outfile = os.path.join(PATH_ANNOTATION, '{genome}', '{genome}.dict')
    params:
        picard  = SOFTWARE['picard']['executable'],
        java    = SOFTWARE['java']['executable'],
        threads = config['execution']['rules']['fasta_dict']['threads'],
        mem     = config['execution']['rules']['fasta_dict']['memory'],
        tempdir = TEMPDIR,
        app_name = 'CreateSequenceDictionary'
    log:
        log = os.path.join(PATH_LOG, '{genome}.fasta_dict.log')
    message:
        """
            Fasta dict:
                input  : {input}
                output : {output}
        """
    run:
        tool = java_tool(params.java, params.threads, params.mem, params.tempdir, params.picard, params.app_name)

        command = ' '.join([
        tool,
        'R=' + str(input.infile),
        'O=' + str(output.outfile),
        '2>' + str(log.log)
        ])
        shell(command)

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
        Rscript = PATH_RSCRIPT
    message:
        """
            Changing GTF id:
                input  : {input}
                output : {output}
        """
    run:
        RunRscript(input, output, params, params.script, 'change_gtf_id.R')


# ----------------------------------------------------------------------------- #
rule gtf_to_refflat:
    input:
        dict = os.path.join(PATH_ANNOTATION, '{genome}', '{genome}.dict'),
        gtf  = os.path.join(PATH_ANNOTATION, '{genome}', '{genome}.gene_id.gtf')
    output:
        outfile = os.path.join(PATH_ANNOTATION, '{genome}', '{genome}.refFlat')
    params:
        threads   = config['execution']['rules']['gtf_to_refflat']['threads'],
        mem       = config['execution']['rules']['gtf_to_refflat']['memory'],
        java      = SOFTWARE['java']['executable'] ,
        droptools = SOFTWARE['droptools']['executable'],
        tempdir   = TEMPDIR,
        app_name  = 'ConvertToRefFlat'
    log:
        log = os.path.join(PATH_LOG, '{genome}.gtf_to_refflat.log')
    message:"""
            GTF To refFlat:
                input
                    dict : {input.dict}
                    gtf  : {input.gtf}
                output : {output}
        """
    run:
        tool = java_tool(params.java, params.threads, params.mem, params.tempdir, params.droptools, params.app_name)

        command = ' '.join([
        tool,
        'ANNOTATIONS_FILE='     + str(input.gtf),
        'SEQUENCE_DICTIONARY='  + str(input.dict),
        'O=' + str(output.outfile),
        '2>' + str(log.log)
        ])
        print(command, file=sys.stderr)
        shell(command)


# # ----------------------------------------------------------------------------- #
def fetch_fastq_for_merge(wc):
    reads = SAMPLE_SHEET.lookup('Merged_BAM', os.path.join(PATH_MAPPED, wc.name, "Fastq_" + wc.num + ".bam"), ['reads'])
    barcode = SAMPLE_SHEET.lookup('Merged_BAM', os.path.join(PATH_MAPPED, wc.name, "Fastq_" + wc.num + ".bam"), ['barcode'])
    h = { 'reads'   : [os.path.join(PATH_FASTQ, file) for file in reads],
          'barcode' : [os.path.join(PATH_FASTQ, file) for file in barcode]}
    return(h)

rule merge_fastq_to_bam:
    input:
        unpack(fetch_fastq_for_merge)
    output:
        outfile = os.path.join(PATH_MAPPED, "{name}", "Fastq_" + "{num}" + ".bam")
    params:
        name    = '{name}' + '_' + '{num}',
        picard  = SOFTWARE['picard']['executable'],
        java    = SOFTWARE['java']['executable'],
        threads = config['execution']['rules']['merge_fastq_to_bam']['threads'],
        mem     = config['execution']['rules']['merge_fastq_to_bam']['memory'],
        tempdir = TEMPDIR,
        app_name = 'FastqToSam'
    log:
        log = os.path.join(PATH_LOG, '{name}.{num}.merge_fastq_to_bam.log')
    message:"""
            merge_fastq_to_bam:
                input:
                    barcode : {input.barcode}
                    reads   : {input.reads}
                output : {output}
        """
    run:
        tool = java_tool(params.java, params.threads, params.mem, params.tempdir, params.picard, params.app_name)

        command = ' '.join([
        tool,
        'O='  + str(output.outfile),
        'F1=' + str(input.barcode),
        'F2=' + str(input.reads),
        'QUALITY_FORMAT=Standard',
        'SAMPLE_NAME=' + str(params.name),
        'SORT_ORDER=unsorted',
        '2>' + str(log.log)
        ])
        print(command, file=sys.stderr)
        shell(command)

# ----------------------------------------------------------------------------- #
def fetch_all_bams(wc):
    bam_files = SAMPLE_SHEET.lookup('sample_name', wc.name, ['Merged_BAM'])
    return(bam_files)

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

# ----------------------------------------------------------------------------- #
rule tag_cells:
    input:
        infile = rules.merge_bam_per_sample.output.outfile
    output:
        outfile   = temp(os.path.join(PATH_MAPPED, "{name}", "{genome}","unaligned_tagged_Cell.bam"))
    params:
        summary   = os.path.join(PATH_MAPPED, "{name}", "{genome}","unaligned_tagged_Cellular.bam_summary.txt"),
        droptools  = SOFTWARE['droptools']['executable'],
        java       = SOFTWARE['java']['executable'],
        app_name   = 'TagBamWithReadSequenceExtended',
        name       = '{name}',
        threads   = config['execution']['rules']['tag_cells']['threads'],
        mem       = config['execution']['rules']['tag_cells']['memory'],
        tempdir    = TEMPDIR,
        base_qual  = 10,
        barcoded_read = 1,
        discard_read  = 'false'
    log:
       log = os.path.join(PATH_LOG, "{name}.{genome}.tag_cells.log")
    message:"""
            tag cells
                input: {input.infile}
                output : {output.outfile}
        """
    run:
        tool = java_tool(params.java, params.threads, params.mem, params.tempdir, params.droptools, params.app_name)

        # fetches method from the sample_sheet
        method = SAMPLE_SHEET.lookup('sample_name', params.name, ['method'])[0]
        adapter_params = ADAPTER_PARAMETERS[method]['cell_barcode']

        command = ' '.join([
        tool,
        'SUMMARY='       + str(params.summary),
        'BASE_RANGE='    + "-".join([str(adapter_params['base_min']),str(adapter_params['base_max'])]),
        'BASE_QUALITY='  + str(params.base_qual),
        'BARCODED_READ=' + str(params.barcoded_read),
        'DISCARD_READ='  + str(params.discard_read),
        'TAG_NAME=XC',
        'NUM_BASES_BELOW_QUALITY=1',
        'INPUT='  + str(input.infile),
        'OUTPUT=' + str(output.outfile)
        ])
        print(command, file=sys.stderr)
        shell(command)

# ----------------------------------------------------------------------------- #
rule tag_molecules:
    input:
        infile = rules.tag_cells.output.outfile
    output:
        outfile   = temp(os.path.join(PATH_MAPPED, "{name}", "{genome}","unaligned_tagged_CellMolecular.bam"))
    params:
        summary   = os.path.join(PATH_MAPPED, "{name}", "{genome}","unaligned_tagged_Molecular.bam_summary.txt"),
        droptools = SOFTWARE['droptools']['executable'],
        java      = SOFTWARE['java']['executable'],
        app_name  = 'TagBamWithReadSequenceExtended',
        name      = '{name}',
        threads   = config['execution']['rules']['tag_molecules']['threads'],
        mem       = config['execution']['rules']['tag_molecules']['memory'],
        tempdir   = TEMPDIR,
        base_qual = 10,
        barcoded_read = 1,
        discard_read  = 'true'
    log:
       log = os.path.join(PATH_LOG, "{name}.{genome}.tag_molecules.log")

    run:
        tool = java_tool(params.java, params.threads, params.mem, params.tempdir, params.droptools, params.app_name)

        method = SAMPLE_SHEET.lookup('sample_name', params.name, ['method'])[0]
        adapter_params = ADAPTER_PARAMETERS[method]['umi_barcode']

        command = ' '.join([
        tool,
        'SUMMARY='       + params.summary,
        'BASE_RANGE='    + "-".join([str(adapter_params['base_min']),str(adapter_params['base_max'])]),
        'BASE_QUALITY='  + str(params.base_qual),
        'BARCODED_READ=' + str(params.barcoded_read),
        'DISCARD_READ='  + str(params.discard_read),
        'TAG_NAME=XM',
        'NUM_BASES_BELOW_QUALITY=1',
        'INPUT='  + str(input.infile),
        'OUTPUT=' + str(output.outfile)
        ])
        shell(command)

# ----------------------------------------------------------------------------- #
rule filter_bam:
    input:
        infile = rules.tag_molecules.output.outfile
    output:
        outfile   = temp(os.path.join(PATH_MAPPED, "{name}", "{genome}","unaligned_tagged_filtered.bam"))
    params:
        droptools = SOFTWARE['droptools']['executable'],
        java      = SOFTWARE['java']['executable'],
        app_name  = 'FilterBAM',
        threads   = config['execution']['rules']['filter_bam']['threads'],
        mem       = config['execution']['rules']['filter_bam']['memory'],
        tempdir   = TEMPDIR
    log:
       log = os.path.join(PATH_LOG, "{name}.{genome}.filter_bam.log")

    run:
        tool = java_tool(params.java, params.threads, params.mem, params.tempdir, params.droptools, params.app_name)

        command = ' '.join([
        tool,
        'TAG_REJECT=XQ',
        'INPUT='  + str(input.infile),
        'OUTPUT=' + str(output.outfile)
        ])
        shell(command)

# ----------------------------------------------------------------------------- #
rule trim_starting_sequence:
    input:
        infile  = rules.filter_bam.output.outfile
    output:
        outfile = temp(os.path.join(PATH_MAPPED, "{name}", "{genome}","unaligned_tagged_trimmed_smart.bam"))
    params:
        droptools  = SOFTWARE['droptools']['executable'],
        java       = SOFTWARE['java']['executable'],
        app_name   = 'TrimStartingSequence',
        threads    = config['execution']['rules']['trim_starting_sequence']['threads'],
        mem        = config['execution']['rules']['trim_starting_sequence']['memory'],
        tempdir    = TEMPDIR,
        summary    = os.path.join(PATH_MAPPED, "{name}", "{genome}","adapter_trimming_report.txt"),
        mismatches = 1,
        num_bases  = 5,
        sequence   = 'AAGCAGTGGTATCAACGCAGAGTGAATGGG'
    log:
       log = os.path.join(PATH_LOG, "{name}.{genome}.trim_starting_sequence.log")

    run:
        tool = java_tool(params.java, params.threads, params.mem, params.tempdir, params.droptools, params.app_name)

        command = ' '.join([
        tool,
        'OUTPUT_SUMMARY='  + str(params.summary),
        'MISMATCHES='      + str(params.mismatches),
        'NUM_BASES='       + str(params.num_bases),
        'SEQUENCE='        + str(params.sequence),
        'INPUT='           + str(input.infile),
        'OUTPUT='          + str(output.outfile)
        ])
        shell(command)

# ----------------------------------------------------------------------------- #
rule trim_polya:
    input:
        infile = rules.trim_starting_sequence.output.outfile
    output:
        outfile   = temp(os.path.join(PATH_MAPPED, "{name}", "{genome}","unaligned_mc_tagged_polyA_filtered.bam"))
    params:
        droptools  = SOFTWARE['droptools']['executable'],
        java       = SOFTWARE['java']['executable'],
        app_name   = 'PolyATrimmer',
        threads    = config['execution']['rules']['trim_polya']['threads'],
        mem        = config['execution']['rules']['trim_polya']['memory'],
        tempdir    = TEMPDIR,
        summary    = os.path.join(PATH_MAPPED, "{name}", "{genome}","polyA_trimming_report.txt"),
        mismatches = 0,
        num_bases  = 6
    log:
       log = os.path.join(PATH_LOG, "{name}.{genome}.trim_polya.log")

    run:
        tool = java_tool(params.java, params.threads, params.mem, params.tempdir, params.droptools, params.app_name)

        command = ' '.join([
        tool,
        'OUTPUT_SUMMARY='  + str(params.summary),
        'MISMATCHES='      + str(params.mismatches),
        'NUM_BASES='       + str(params.num_bases),
        'INPUT='           + str(input.infile),
        'OUTPUT='          + str(output.outfile)
        ])
        shell(command)
# ----------------------------------------------------------------------------- #
rule sam_to_fastq:
    input:
        infile = rules.trim_polya.output.outfile
    output:
        outfile   = temp(os.path.join(PATH_MAPPED, "{name}", "{genome}","unaligned_mc_tagged_polyA_filtered.fastq"))
    params:
        picard     = SOFTWARE['picard']['executable'],
        java       = SOFTWARE['java']['executable'],
        app_name   = 'SamToFastq',
        threads    = config['execution']['rules']['sam_to_fastq']['threads'],
        mem        = config['execution']['rules']['sam_to_fastq']['memory'],
        tempdir    = TEMPDIR
    log:
       log = os.path.join(PATH_LOG, "{name}.{genome}.SamToFastq.log")

    run:
        tool = java_tool(params.java, params.threads, params.mem, params.tempdir, params.picard, params.app_name)

        command = ' '.join([
        tool,
        'INPUT='   + str(input.infile),
        'FASTQ='   + str(output.outfile)
        ])
        shell(command)

# ----------------------------------------------------------------------------- #
rule map_star:
    input:
        infile = rules.sam_to_fastq.output.outfile,
        genome = rules.make_star_reference.output
    output:
        outfile   = temp(os.path.join(PATH_MAPPED, "{name}", "{genome}","star.Aligned.out.sam"))
    params:
        star       = SOFTWARE['star']['executable'],
        star_args  = SOFTWARE['star']['args'],
        genome     = os.path.join(PATH_ANNOTATION, '{genome}','STAR_INDEX'),
        outpath    = os.path.join(PATH_MAPPED, "{name}", "{genome}"),
        threads    = config['execution']['rules']['map_star']['threads'],
        mem        = config['execution']['rules']['map_star']['memory'],
        tempdir    = TEMPDIR,
    log:
       log = os.path.join(PATH_LOG, "{name}.{genome}.star.log")

    shell:"""
    {params.star} {params.star_args} --genomeDir {params.genome} --runThreadN {params.threads} --outFileNamePrefix {params.outpath}/star. --readFilesIn {input.infile}
    """

# ----------------------------------------------------------------------------- #
rule sort_aligned:
    input:
        infile = rules.map_star.output.outfile
    output:
        outfile   = temp(os.path.join(PATH_MAPPED, "{name}", "{genome}","aligned.sorted.bam"))
    params:
        picard     = SOFTWARE['picard']['executable'],
        java       = SOFTWARE['java']['executable'],
        app_name   = 'SortSam',
        threads    = config['execution']['rules']['sort_aligned']['threads'],
        mem        = config['execution']['rules']['sort_aligned']['memory'],
        tempdir    = TEMPDIR
    log:
       log = os.path.join(PATH_LOG, "{name}.{genome}.sort_aligned.log")

    run:
        tool = java_tool(params.java, params.threads, params.mem, params.tempdir, params.picard, params.app_name)

        command = ' '.join([
        tool,
        'SORT_ORDER=queryname',
        'TMP_DIR=' + str(params.tempdir),
        'INPUT='   + str(input.infile),
        'OUTPUT='  + str(output.outfile)
        ])
        shell(command)

# ----------------------------------------------------------------------------- #
rule merge_bam:
    input:
        mapped    = rules.sort_aligned.output.outfile,
        unmapped  = rules.trim_polya.output.outfile,
        reference = os.path.join(PATH_ANNOTATION, '{genome}', '{genome}.fasta'),
        dict      = rules.fasta_dict.output.outfile
    output:
        outfile   = temp(os.path.join(PATH_MAPPED, "{name}", "{genome}","merged.bam"))
    params:
        picard     = SOFTWARE['picard']['executable'],
        java       = SOFTWARE['java']['executable'],
        app_name   = 'MergeBamAlignment',
        threads    = config['execution']['rules']['merge_bam']['threads'],
        mem        = config['execution']['rules']['merge_bam']['memory'],
        tempdir    = TEMPDIR
    log:
       log = os.path.join(PATH_LOG, "{name}.{genome}.merge_bam.log")

    run:
        tool = java_tool(params.java, params.threads, params.mem, params.tempdir, params.picard, params.app_name)

        command = ' '.join([
        tool,
        'REFERENCE_SEQUENCE=' + str(input.reference),
        'UNMAPPED_BAM='       + str(input.unmapped),
        'ALIGNED_BAM='        + str(input.mapped),
        'INCLUDE_SECONDARY_ALIGNMENTS=false',
        'PAIRED_RUN=false',
        'OUTPUT='             + str(output.outfile)
        ])
        shell(command)

# ----------------------------------------------------------------------------- #
rule tag_with_gene_exon:
    input:
        infile    = rules.merge_bam.output.outfile,
        refflat   = rules.gtf_to_refflat.output
    output:
        outfile   = os.path.join(PATH_MAPPED, "{name}", "{genome}","star_gene_exon_tagged.bam")
    params:
        droptools  = SOFTWARE['droptools']['executable'],
        java       = SOFTWARE['java']['executable'],
        app_name   = 'TagReadWithGeneExon',
        threads    = config['execution']['rules']['tag_with_gene_exon']['threads'],
        mem        = config['execution']['rules']['tag_with_gene_exon']['memory'],
        tempdir    = TEMPDIR
    log:
       log = os.path.join(PATH_LOG, "{name}.{genome}.tag_with_gene_exon.log")

    run:
        tool = java_tool(params.java, params.threads, params.mem, params.tempdir, params.droptools, params.app_name)

        command = ' '.join([
        tool,
        'ANNOTATIONS_FILE=' + str(input.refflat),
        'TAG=GE',
        'CREATE_INDEX=true',
        'INPUT='  + str(input.infile),
        'OUTPUT=' + str(output.outfile)
        ])
        shell(command)



# ----------------------------------------------------------------------------- #
rule extract_read_statistics:
    input:
        bamfile = rules.tag_with_gene_exon.output.outfile
    output:
        outfile = os.path.join(PATH_MAPPED, "{name}", "{genome}",'{name}_{genome}_ReadStatistics.txt')
    params:
        outname  = "{name}_{genome}",
        threads  = config['execution']['rules']['extract_read_statistics']['threads'],
        mem      = config['execution']['rules']['extract_read_statistics']['memory'],
        script   = PATH_SCRIPT,
        Rscript  = PATH_RSCRIPT
    message: """
            extract_read_statistics:
                input:  {input.bamfile}
                output: {output.outfile}
        """
    run:
        RunRscript(input, output, params, params.script, 'Extract_Read_Statistics.R')




# ----------------------------------------------------------------------------- #
# calculates the number of reads per cell
rule bam_tag_histogram:
    input:
        infile = rules.tag_with_gene_exon.output.outfile
    output:
        outfile = os.path.join(PATH_MAPPED, "{name}", "{genome}",'{name}_{genome}_BAMTagHistogram.txt')
    params:
        outdir    = os.path.join(PATH_MAPPED, "{name}", "{genome}"),
        outname   = "{name}_{genome}",
        threads   = config['execution']['rules']['bam_tag_histogram']['threads'],
        mem       = config['execution']['rules']['bam_tag_histogram']['memory'],
        java      = SOFTWARE['java']['executable'],
        droptools = SOFTWARE['droptools']['executable'],
        tempdir   = TEMPDIR,
        app_name  = 'BAMTagHistogram'
    message: """
            BamTagHistogram:
                input:  {input.infile}
                output: {output.outfile}
        """
    run:
        tool = java_tool(params.java, params.threads, params.mem, params.tempdir, params.droptools, params.app_name)

        command = ' '.join([
        tool,
        'O=' + str(output.outfile),
        'I=' + str(input.infile),
        'TAG=XC'
        ])
        shell(command)


# ----------------------------------------------------------------------------- #
# calculates the UMI matrix
rule find_absolute_read_cutoff:
    input:
        infile = rules.bam_tag_histogram.output.outfile
    output:
        outfile = os.path.join(PATH_MAPPED, "{name}", "{genome}",'{name}_{genome}_ReadCutoff.yaml')
    params:
        outdir   = os.path.join(PATH_MAPPED, "{name}", "{genome}"),
        outname  = "{name}_{genome}",
        threads  = config['execution']['rules']['find_absolute_read_cutoff']['threads'],
        mem      = config['execution']['rules']['find_absolute_read_cutoff']['memory'],
        cutoff   = 50000,
        script   = PATH_SCRIPT,
        Rscript  = PATH_RSCRIPT
    message: """
            find_absolute_read_cutoff:
                input:  {input.infile}
                output: {output.outfile}
        """
    run:
        RunRscript(input, output, params, params.script, 'Find_Absolute_Read_Cutoff.R')


# ----------------------------------------------------------------------------- #
# calculates the UMI matrix
rule get_umi_matrix:
    input:
        infile        = rules.tag_with_gene_exon.output,
        reads_cutoff  = rules.find_absolute_read_cutoff.output.outfile
    output:
        outfile = os.path.join(PATH_MAPPED, "{name}", "{genome}",'{name}_{genome}_UMI.Matrix.txt')
    params:
        outdir            = os.path.join(PATH_MAPPED, "{name}", "{genome}"),
        outname           = "{name}_{genome}",
        threads           = config['execution']['rules']['get_umi_matrix']['threads'],
        mem               = config['execution']['rules']['get_umi_matrix']['memory'],
        java              = SOFTWARE['java']['executable'] ,
        droptools         = SOFTWARE['droptools']['executable'],
        app_name          = 'DigitalExpression',
        tempdir           = TEMPDIR,
        genes_per_cell    = SOFTWARE['droptools']['genes_per_cell'],
        num_core_barcodes = SOFTWARE['droptools']['num_core_barcodes']
    message: """
            Count UMI:
                input:  {input.infile}
                reads:  {input.reads_cutoff}
                output: {output.outfile}
        """
    run:
        with open(input.reads_cutoff) as stream:
            reads_cutoff = yaml.load(stream)['reads_cutoff']

        tool = java_tool(params.java, params.threads, params.mem, params.tempdir, params.droptools, params.app_name)

        command = ' '.join([
        tool,
        'O=' + str(output.outfile),
        'I=' + str(input.infile),
        'SUMMARY=' + os.path.join(params.outdir, params.outname + '_Summary.txt'),
#         'MIN_NUM_GENES_PER_CELL=' + str(params.genes_per_cell),
        'MIN_NUM_READS_PER_CELL=' + str(reads_cutoff),
        'OUTPUT_READS_INSTEAD=F'
        ])
        shell(command)


# ----------------------------------------------------------------------------- #
# calculates the reads matrix - used for PCR duplication estimations
rule get_reads_matrix:
    input:
        infile        = rules.tag_with_gene_exon.output,
        reads_cutoff  = rules.find_absolute_read_cutoff.output.outfile
    output:
        outfile = os.path.join(PATH_MAPPED, "{name}", "{genome}",'{name}_{genome}_READS.Matrix.txt')
    params:
        outdir            = os.path.join(PATH_MAPPED, "{name}", "{genome}"),
        outname           = "{name}_{genome}",
        java              = SOFTWARE['java']['executable'] ,
        droptools         = SOFTWARE['droptools']['executable'],
        app_name          = 'DigitalExpression',
        tempdir           = TEMPDIR,
        threads           = config['execution']['rules']['get_reads_matrix']['threads'],
        mem               = config['execution']['rules']['get_reads_matrix']['memory']
    message: """
            Count UMI:
                input:  {input.infile}
                reads:  {input.reads_cutoff}
                output: {output.outfile}
        """
    run:
        with open(input.reads_cutoff) as stream:
            reads_cutoff = yaml.load(stream)['reads_cutoff']

        tool = java_tool(params.java, params.threads, params.mem, params.tempdir, params.droptools, params.app_name)

        command = ' '.join([
        tool,
        'O=' + str(output.outfile),
        'I=' + str(input.infile),
        'SUMMARY=' + os.path.join(params.outdir, params.outname + '_Summary.txt'),
        'MIN_NUM_READS_PER_CELL=' + str(reads_cutoff),
        'OUTPUT_READS_INSTEAD=T'
        ])
        shell(command)

# ----------------------------------------------------------------------------- #
rule extract_downstream_statistics:
    input:
        umi_matrix   = rules.get_umi_matrix.output.outfile,
        reads_matrix = rules.get_reads_matrix.output.outfile,
        reads_stats  = rules.extract_read_statistics.output.outfile,
        reads_cutoff = rules.find_absolute_read_cutoff.output.outfile
    output:
        outfile = os.path.join(PATH_MAPPED, "{name}", "{genome}",'{name}_{genome}_DownstreamStatistics.txt')
    params:
        file_location = os.path.join(PATH_MAPPED, "{name}", "{genome}"),
        outname       = "{name}_{genome}",
        threads       = config['execution']['rules']['extract_downstream_statistics']['threads'],
        mem           = config['execution']['rules']['extract_downstream_statistics']['memory'],
        script        = PATH_SCRIPT,
        Rscript       = PATH_RSCRIPT
    message: """
            extract_downstream_statistics:
                umi:    {input.umi_matrix}
                reads:  {input.reads_matrix}
                stats:  {input.reads_stats}
                output: {output.outfile}
        """
    run:
        RunRscript(input, output, params, params.script, 'Extract_Downstream_Statistics.R')

# ----------------------------------------------------------------------------- #
# convert UMI matrix from txt format into one loom format
rule convert_matrix_from_txt_to_loom:
    input:
        infile        = os.path.join(PATH_MAPPED, "{name}", "{genome}",'{name}_{genome}_UMI.Matrix.txt'),
        gtf           = lambda wildcards: os.path.join(PATH_ANNOTATION, wildcards.genome, '.'.join([wildcards.genome, 'gtf']))
    output:
        outfile       = os.path.join(PATH_MAPPED, "{name}", "{genome}",'{name}_{genome}_UMI.Matrix.loom')
    params:
        python      = SOFTWARE['python']['executable'],
        python_args = SOFTWARE['python']['args']
    log:
        log = os.path.join(PATH_LOG, "{name}.{genome}.convert2loom.log")
    message: """
            Convert UMI Matrix from .txt to .loom:
                input:  {input.infile}
                output: {output.outfile}
        """
    shell: "{params.python} {params.python_args} {PATH_SCRIPT}/convert_matrix_to_loom.py {wildcards.name} {input.infile} {input.gtf} {output.outfile} &> {log.log}"


# ----------------------------------------------------------------------------- #
## combines multiple loom files into one loom file
rule combine_UMI_matrices_into_loom:
    input:
        infile        = lambda wildcards: expand(os.path.join(PATH_MAPPED, "{name}", wildcards.genome, '_'.join(["{name}", wildcards.genome, 'UMI.Matrix.loom'])), name = SAMPLE_NAMES)
    output:
        outfile       = os.path.join(PATH_MAPPED, "{genome}_UMI.loom")
    params:
         tmpfile = os.path.join(PATH_MAPPED, "{genome}_UMI.loom.tmp"),
         python  = SOFTWARE['python']['executable'],
         python_args = SOFTWARE['python']['args'],
         rm      = SOFTWARE['rm']['executable']
    log:
        log = os.path.join(PATH_LOG, "{genome}.combine_looms.log")
    message: """
            Combine multiple loom files into one:
                input:  {input.infile}
                output: {output.outfile}
        """
    shell: "echo {input.infile} > {params.tmpfile}; {params.python} {params.python_args} {PATH_SCRIPT}/combine_UMI_matrices.py {params.tmpfile} {output.outfile}; {params.rm} {params.tmpfile} &> {log.log}"


# ----------------------------------------------------------------------------- #
## Imports and preprocesses the combined loom files and saves as SingleCellExperiment.RDS objects.
rule convert_loom_to_singleCellExperiment:
    input:
        infile        = os.path.join(PATH_MAPPED, "{genome}_UMI.loom")
    output:
        outfile       = os.path.join(PATH_MAPPED, "{genome}.SingleCellExperiment.RDS")
    log:
        log = os.path.join(PATH_LOG, "{genome}.loom2sce.log")
    params:
        Rscript = PATH_RSCRIPT
    message: """
            Import loom file, preprocess and save as singleCellExperiment.RDS:
                input:  {input.infile}
                output: {output.outfile}
        """
    shell: "{params.Rscript} --vanilla {PATH_SCRIPT}/convert_loom_to_singleCellExperiment.R --loomFile={input.infile} --sampleSheetFile={PATH_SAMPLE_SHEET} --gtfFile={PATH_GTF_PRIMARY} --genomeBuild={wildcards.genome} --outFile={output.outfile} &> {log.log}"


# ----------------------------------------------------------------------------- #
## Using the preprocessed SingleCellExperiment.RDS file, generates a self-contained HTML report
rule report:
    input:
        infile        = os.path.abspath(os.path.join(PATH_MAPPED, "{genome}.SingleCellExperiment.RDS")),
        read_stats    = READ_STATISTICS
    output:
        outfile       = os.path.join(PATH_MAPPED, "{genome}.scRNA-Seq.report.html")
    params:
        reportRmd     = os.path.join(PATH_SCRIPT, "scrnaReport.Rmd"),
        Rscript       = PATH_RSCRIPT,
        workdir       = os.path.abspath(PATH_MAPPED),
        sceRdsFile    = lambda wildcards: os.path.abspath(os.path.join(PATH_MAPPED, wildcards.genome + ".SingleCellExperiment.RDS")),
        path_mapped   = PATH_MAPPED
    log:
        log = os.path.join(PATH_LOG, "{genome}.scRNA-Seq.report.log")
    message: """
            Generate an HTML report from SingleCellExperiment.RDS:
                input:  {input.infile}
                output: {output.outfile}
        """
    shell: "{params.Rscript} --vanilla {PATH_SCRIPT}/renderReport.R --reportFile={params.reportRmd} --sceRdsFile={params.sceRdsFile} --covariates='{COVARIATES}' --prefix={wildcards.genome} --workdir={params.workdir} --path_mapped='{params.path_mapped}' 2>&1 {log.log}"

# ----------------------------------------------------------------------------- #
rule bam_to_BigWig:
    input:
        bamfile            = rules.tag_with_gene_exon.output,
        cell_cutoff_file   = rules.find_absolute_read_cutoff.output.outfile,
        reads_by_cell_file = rules.bam_tag_histogram.output.outfile
    output:
        bwfile = os.path.join(PATH_MAPPED, "{name}", "{genome}",'{name}_{genome}.bw')
    params:
        threads = config['execution']['rules']['bam_to_BigWig']['threads'],
        mem     = config['execution']['rules']['bam_to_BigWig']['memory'],
        script  = PATH_SCRIPT,
        Rscript = PATH_RSCRIPT
    message: """
            bam_to_BigWig:
                input:  {input.bamfile}
                output: {output.bwfile}
            """
    run:
        RunRscript(input, output, params, params.script, 'BamToBigWig.R')


# ----------------------------------------------------------------------------- #
rule fastqc:
    input:
        infile = os.path.join(PATH_FASTQ, "{name}")
    output:
        outfile =  os.path.join(PATH_FASTQC, "{name}.fastqc.done")
    params:
        outpath = os.path.join(PATH_FASTQC),
        threads = config['execution']['rules']['fastqc']['threads'],
        mem     = config['execution']['rules']['fastqc']['memory'],
        java    = SOFTWARE['java']['executable'],
        fastqc  = SOFTWARE['fastqc']['executable']
    log:
        log = os.path.join(PATH_LOG, "{name}.fastqc.log")
    message: """
            fastqc:
                input: {input.infile}
                output: {output.outfile}
            """
    shell:"""
        {params.fastqc} -j {params.java} -t {params.threads} -o {params.outpath} {input.infile} 2> {log.log}
        touch {output.outfile}
    """
