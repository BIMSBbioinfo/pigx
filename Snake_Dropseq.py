
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

PATH_SCRIPT = os.path.join(config['locations']['pkglibexecdir'])

# loads function for running Rscripts
include: os.path.join(PATH_SCRIPT, 'Run_Rscript.py')
# functions for java input
include: os.path.join(PATH_SCRIPT, 'Accessory_Functions.py')
include: os.path.join(PATH_SCRIPT, 'validate_input.py')
validate_config(config)


# ----------------------------------------------------------------------------- #
# Software parameters
SOFTWARE = config['tools']
# variables

GENOME_NAME_PRIMARY = config['annotation']['primary']['genome']['name']
REFERENCE_NAMES = [GENOME_NAME_PRIMARY]
COVARIATES = config['covariates']

# ----------------------------------------------------------------------------- #
# PATHS
OUTPUT_DIR = config['locations']['output-dir']
PATH_FASTQ = config['locations']['reads-dir']
TEMPDIR    = config['locations']['tempdir']
PATH_ANNOTATION = os.path.join(OUTPUT_DIR, 'Annotation')
PATH_MAPPED   = os.path.join(OUTPUT_DIR, 'Mapped')
PATH_LOG      = os.path.join(OUTPUT_DIR, 'Log')
SAMPLE_SHEET_FILE = config['locations']['sample-sheet']
PATH_RSCRIPT = SOFTWARE['Rscript']['executable']

PATH_ANNOTATION_PRIMARY = os.path.join(PATH_ANNOTATION, GENOME_NAME_PRIMARY)
PATH_REFERENCE_PRIMARY  = config['annotation']['primary']['genome']['fasta']
PATH_GTF_PRIMARY        = config['annotation']['primary']['gtf']
PATH_META_DATA          = config['locations']['metadata']


GENOME_NAME_MIX = None
if 'secondary' in set(config['annotation'].keys()):
    GENOME_NAME_SECONDARY    = config['annotation']['secondary']['genome']['name']
    PATH_REFERENCE_SECONDARY = config['annotation']['secondary']['genome']['fasta']
    PATH_GTF_SECONDARY       = config['annotation']['secondary']['gtf']

    GENOME_NAME_MIX = GENOME_NAME_PRIMARY + '_' + GENOME_NAME_SECONDARY
    PATH_ANNOTATION_MIX = os.path.join(PATH_ANNOTATION, GENOME_NAME_MIX)
    REFERENCE_NAMES = REFERENCE_NAMES + [GENOME_NAME_MIX]

## Load sample sheet
with open(SAMPLE_SHEET_FILE, 'r') as fp:
  rows =  [row for row in csv.reader(fp, delimiter=',')]
  header = rows[0]; rows = rows[1:]
  SAMPLE_SHEET = [dict(zip(header, row)) for row in rows]

# Convenience function to access fields of sample sheet columns that
# match the predicate.  The predicate may be a string.
def lookup(column, predicate, fields=[]):
  if inspect.isfunction(predicate):
    records = [line for line in SAMPLE_SHEET if predicate(line[column])]
  else:
    records = [line for line in SAMPLE_SHEET if line[column]==predicate]
  return [record[field] for record in records for field in fields]

SAMPLE_NAMES = [line['sample_name'] for line in SAMPLE_SHEET]

# ----------------------------------------------------------------------------- #
# sets the temporrary directory to default in the working directory if the tempdir does not exist
if TEMPDIR == None:
    if TMPDIR in os.environ.keys():
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
FASTQC = expand(os.path.join(PATH_MAPPED, "{name}", "{name}.fastqc.done"), name=SAMPLE_NAMES)

# ----------------------------------------------------------------------------- #
# Change reference gene_name to gene_id
PATH_GTF_PRIMARY_ID = expand(os.path.join(PATH_ANNOTATION_PRIMARY, "{name}" + "gene_id.gtf"), name=REFERENCE_NAMES)

# ----------------------------------------------------------------------------- #
# STAR INDEX
MAKE_STAR_INDEX = expand(os.path.join(PATH_ANNOTATION, '{genome}','STAR_INDEX','done.txt'), genome = REFERENCE_NAMES)

# ----------------------------------------------------------------------------- #
# MERGE BARODE AND READS FASTQ FILES
MERGE_FASTQ_TO_BAM = expand(os.path.join(PATH_MAPPED, "{name}", "{name}" + '.fastq.bam'), name=SAMPLE_NAMES)

# ----------------------------------------------------------------------------- #
# MAPPING
MAP_scRNA = expand(os.path.join(PATH_MAPPED, "{name}", "{genome}", "star_gene_exon_tagged.bam"), genome = REFERENCE_NAMES, name = SAMPLE_NAMES)

# ----------------------------------------------------------------------------- #
# READ statistics
READ_STATISTICS = expand(os.path.join(PATH_MAPPED, "{name}", "{genome}",'{name}_{genome}_ReadStatistics.txt'), genome = REFERENCE_NAMES, name = SAMPLE_NAMES)
# ----------------------------------------------------------------------------- #
# Number of reads per cell calculation
BAM_HISTOGRAM = expand(os.path.join(PATH_MAPPED, "{name}", "{genome}",'{name}_{genome}_BAMTagHistogram.txt'), genome = REFERENCE_NAMES, name = SAMPLE_NAMES)

# ----------------------------------------------------------------------------- #
# Number of reads per cell calculation
FIND_READ_CUTOFF = expand(os.path.join(PATH_MAPPED, "{name}", "{genome}",'{name}_{genome}_ReadCutoff.yaml'), genome = REFERENCE_NAMES, name = SAMPLE_NAMES)

# ----------------------------------------------------------------------------- #
# UMI matrix
UMI = expand(os.path.join(PATH_MAPPED, "{name}", "{genome}",'{name}_{genome}_UMI.Matrix.txt'), genome = REFERENCE_NAMES, name = SAMPLE_NAMES)

# UMI matrix in loom format
UMI_LOOM =  expand(os.path.join(PATH_MAPPED, "{name}", "{genome}",'{name}_{genome}_UMI.Matrix.loom'), genome = REFERENCE_NAMES, name = SAMPLE_NAMES)

# Combined UMI matrices in loom format
COMBINED_UMI_MATRICES = expand(os.path.join(PATH_MAPPED, "{genome}_UMI.loom"), genome = REFERENCE_NAMES)

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

RULE_ALL = RULE_ALL + DICT + REFFLAT + MAKE_STAR_INDEX + FASTQC + MERGE_FASTQ_TO_BAM + MAP_scRNA + READ_STATISTICS + BAM_HISTOGRAM + FIND_READ_CUTOFF + UMI + BIGWIG + UMI_LOOM + COMBINED_UMI_MATRICES + SCE_RDS_FILES + REPORT_FILES

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
        fasta = LINK_REFERENCE_PRIMARY,
    params:
        threads=1,
        mem='4G'
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
            threads=1,
            mem='4G',
            genome_name_primary   = GENOME_NAME_PRIMARY,
            genome_name_secondary = GENOME_NAME_SECONDARY
        message:
            """
                Combining fasta files:
                    primary   : {input.primary}
                    secondary : {input.secondary}
                    output: {output.outfile}
            """
        shell:"""
            cat {input.primary}   | perl -pe 's|^>|>{params.genome_name_primary}|' >     {output.outfile}
            cat {input.secondary} | perl -pe 's|^>|>{params.genome_name_secondary}|' >> {output.outfile}
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
        star    = SOFTWARE['star'],
        threads = 8,
        mem     = '40G'
    log:
        os.path.join(PATH_LOG, '{genome}.make_star_reference.log')
    message:"""
        Star reference:
            input:
                fasta : {input.fasta}
                gtf   : {input.gtf}
        """
    shell:"""
        {params.star} --runMode genomeGenerate --genomeDir {params.outdir} --genomeFastaFiles {input.fasta} --runThreadN {params.threads} --sjdbGTFfile {input.gtf} --sjdbOverhang 99
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
            threads=1,
            mem='4G',
            genome_name_primary   = GENOME_NAME_PRIMARY,
            genome_name_secondary = GENOME_NAME_SECONDARY
        message:
            """
                Combining gtf files:
                    primary   : {input.primary}
                    secondary : {input.secondary}
                    output: {output.outfile}
            """
        shell:"""
            cat {input.primary}   | perl -pe 's|^|{params.genome_name_primary}|' > {output.outfile}
            cat {input.secondary} | perl -pe 's|^|{params.genome_name_secondary}|' >> {output.outfile}
    """

# ----------------------------------------------------------------------------- #

rule fasta_dict:
    input:
        os.path.join(PATH_ANNOTATION, '{genome}', '{genome}.fasta')
    output:
        os.path.join(PATH_ANNOTATION, '{genome}', '{genome}.dict')
    params:
        picard  = SOFTWARE['picard']['executable'],
        java    = SOFTWARE['java']['executable']
        threads = 1,
        mem     = '2G',
        tempdir = TEMPDIR
    log:
        os.path.join(PATH_LOG, '{genome}.fasta_dict.log')
    message:
        """
            Fasta dict:
                input  : {input}
                output : {output}
        """
    shell:"""
        {params.java} -XX:ParallelGCThreads={params.threads} -Xmx${params.mem} -Djava.io.tmpdir={params.tempdir} -jar {params.picard} CreateSequenceDictionary R={input} O={output} 2> {log}
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
        threads = 1,
        mem     = '4G',
        Rscript = PATH_RSCRIPT
    message:
        """
            Changing GTF id:
                input  : {input}
                output : {output}
        """
    run:
        RunRscript(input, output, params, 'change_gtf_id.R')




    # shell:"""
    #     cat {input} | perl -pe '/gene_id "([A-Z0-9]+?)";/; $gene_id = $1; s/gene_name ".+?"/gene_name "$gene_id"/;' > {output}
    # """

# ----------------------------------------------------------------------------- #
rule gtf_to_refflat:
    input:
        dict = os.path.join(PATH_ANNOTATION, '{genome}', '{genome}.dict'),
        gtf  = os.path.join(PATH_ANNOTATION, '{genome}', '{genome}.gene_id.gtf')
    output:
        os.path.join(PATH_ANNOTATION, '{genome}', '{genome}.refFlat')
    params:
        threads   = 1,
        mem       = '50G'
        java      = SOFTWARE['java']['executable'] ,
        droptools = SOFTWARE['droptools']['executable'],
        tempdir   = TEMPDIR
    log:
        os.path.join(PATH_LOG, '{genome}.gtf_to_refflat.log')
    message:"""
            GTF To refFlat:
                input
                    dict : {input.dict}
                    gtf  : {input.gtf}
                output : {output}
        """
    shell:"""
        {params.java} -XX:ParallelGCThreads={params.threads} -Xmx${params.mem} -Djava.io.tmpdir={params.tempdir} {params.droptools} ConvertToRefFlat  O={output} ANNOTATIONS_FILE={input.gtf} SEQUENCE_DICTIONARY={input.dict} 2> {log} O={output} ANNOTATIONS_FILE={input.gtf} SEQUENCE_DICTIONARY={input.dict} 2> {log}
    """


# # ----------------------------------------------------------------------------- #
def get_fastq_files(wc):
    h = {'barcode' : os.path.join(PATH_FASTQ, lookup('sample_name', wc.name, ['barcode'])[0]),
         'reads'   : os.path.join(PATH_FASTQ, lookup('sample_name', wc.name, ['reads'])[0])}
    return h


rule merge_fastq_to_bam:
    input:
        unpack(get_fastq_files)
    output:
        os.path.join(PATH_MAPPED, "{name}", "{name}.fastq.bam")
    params:
        name    = '{name}',
        picard  = SOFTWARE['picard']['executable'],
        java    = SOFTWARE['picard']['java'],
        threads = 1,
        mem     = '2G',
        tempdir = TEMPDIR
    log:
        os.path.join(PATH_LOG, '{name}.merge_fastq_to_bam.log')
    message:"""
            Merge fastq barcode and reads:
                input:
                    barcode : {input.barcode}
                    reads   : {input.reads}
                output : {output}
        """
    shell: """
    {params.java} -XX:ParallelGCThreads={params.threads} -Xmx${params.mem} -Djava.io.tmpdir={params.tempdir} -jar {params.picard} FastqToSam O={output} F1={input.barcode} F2={input.reads} QUALITY_FORMAT=Standard SAMPLE_NAME={params.name} SORT_ORDER=queryname 2> {log}
    """
# ----------------------------------------------------------------------------- #
rule tag_cells:
    input:
        rules.merge_fastq_to_bam.output.outfile
    output:
        outfile   = temp(os.path.join(PATH_MAPPED, "{name}", "{genome}","unaligned_tagged_Cell.bam"))
    params:
        summary   = os.path.join(PATH_MAPPED, "{name}", "{genome}","unaligned_tagged_Cellular.bam_summary.txt")
        droptools = SOFTWARE['droptools']['executable'],
        java      = SOFTWARE['java']['executable'],
        app_name  = 'TagBamWithReadSequenceExtended',
        threads   = 8,
        mem       = '35G',
        tempdir   = TEMPDIR,
        base_min  = 1,
        base_max  = 12,
        base_qual = 10,
        barcoded_read = 1,
        discard_read  = 'false'
    log:
       log = os.path.join(PATH_LOG, "{name}.{genome}.tag_cells.log")

    run:
        tool = java_tool(params.java, params.threads, params.mem, params.tempdir, params.droptools, params.app_name)


        command = ' '.join([
        tool,
        'SUMMARY='       + str(params.summary),
        'BASE_RANGE='    + "_".join([str(params.base_min),str(params.base_max)]),
        'BASE_QUALITY='  + str(params.base_qual),
        'BARCODED_READ=' + str(params.barcoded_read),
        'DISCARD_READ='  + str(params.discard_read),
        'TAG_NAME=XC',
        'NUM_BASES_BELOW_QUALITY=1',
        'INPUT='  + str(input.infile),
        'OUTPUT=' + str(output.outfile)
        ])
    shell(command)

# ----------------------------------------------------------------------------- #
rule tag_molecules:
    input:
        rules.tag_cells.outfile
    output:
        outfile   = temp(os.path.join(PATH_MAPPED, "{name}", "{genome}","unaligned_tagged_CellMolecular.bam"))
    params:
        summary   = os.path.join(PATH_MAPPED, "{name}", "{genome}","unaligned_tagged_Molecular.bam_summary.txt"),
        droptools = SOFTWARE['droptools']['executable'],
        java      = SOFTWARE['java']['executable'],
        app_name  = 'TagBamWithReadSequenceExtended',
        threads   = 8,
        mem       = '35G',
        tempdir   = TEMPDIR,
        base_min  = 13,
        base_max  = 20,
        base_qual = 10,
        barcoded_read = 1,
        discard_read  = 'true'
    log:
       log = os.path.join(PATH_LOG, "{name}.{genome}.tag_cells.log")

    run:
        tool = java_tool(params.java, params.threads, params.mem, params.tempdir, params.droptools, params.app_name)


        command = ' '.join([
        tool,
        'SUMMARY='       + params.summary,
        'BASE_RANGE='    + "_".join([str(params.base_min),str(params.base_max)]),
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
        rules.tag_molecules.outfile
    output:
        outfile   = temp(os.path.join(PATH_MAPPED, "{name}", "{genome}","unaligned_tagged_filtered.bam"))
    params:
        droptools = SOFTWARE['droptools']['executable'],
        java      = SOFTWARE['java']['executable'],
        app_name  = 'FilterBAM',
        threads   = 1,
        mem       = '10G',
        tempdir   = TEMPDIR,
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
        rules.filter_bam.outfile
    output:
        outfile   = temp(os.path.join(PATH_MAPPED, "{name}", "{genome}","unaligned_tagged_trimmed_smart.bam"))
    params:
        droptools  = SOFTWARE['droptools']['executable'],
        java       = SOFTWARE['java']['executable'],
        app_name   = 'TrimStartingSequence',
        threads    = 1,
        mem        = '10G',
        tempdir    = TEMPDIR,
        summary    = os.path.join(PATH_MAPPED, "{name}", "{genome}","adapter_trimming_report.txt"),
        mismatches = 1,
        num_bases  = 5
    log:
       log = os.path.join(PATH_LOG, "{name}.{genome}.trim_starting_sequence.log")

    run:
        tool = java_tool(params.java, params.threads, params.mem, params.tempdir, params.droptools, params.app_name)

        command = ' '.join([
        tool,
        'OUTPUT_SUMMARY=', + str(param.summary),
        'MISMATCHES=',     + str(params.mismatches),
        'NUM_BASES='       + str(params.num_bases),
        'INPUT='           + str(input.infile),
        'OUTPUT='          + str(output.outfile)
        ])
    shell(command)

# ----------------------------------------------------------------------------- #
rule trim_polya:
    input:
        rules.trim_starting_sequence.outfile
    output:
        outfile   = temp(os.path.join(PATH_MAPPED, "{name}", "{genome}","unaligned_mc_tagged_polyA_filtered.bam"))
    params:
        droptools  = SOFTWARE['droptools']['executable'],
        java       = SOFTWARE['java']['executable'],
        app_name   = 'PolyATrimmer',
        threads    = 1,
        mem        = '10G',
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
        'OUTPUT_SUMMARY=', + str(param.summary),
        'MISMATCHES=',     + str(params.mismatches),
        'NUM_BASES='       + str(params.num_bases),
        'INPUT='           + str(input.infile),
        'OUTPUT='          + str(output.outfile)
        ])
    shell(command)
# ----------------------------------------------------------------------------- #
rule sam_to_fastq:
    input:
        rules.trim_polya.outfile
    output:
        outfile   = temp(os.path.join(PATH_MAPPED, "{name}", "{genome}","unaligned_mc_tagged_polyA_filtered.fastq"))
    params:
        picard     = SOFTWARE['picard']['executable'],
        java       = SOFTWARE['java']['executable'],
        app_name   = 'SamToFastq',
        threads    = 1,
        mem        = '10G',
        tempdir    = TEMPDIR
    log:
       log = os.path.join(PATH_LOG, "{name}.{genome}.SamToFastq.log")

    run:
        tool = java_tool(params.java, params.threads, params.mem, params.tempdir, params.picard, params.app_name)

        command = ' '.join([
        tool,
        'INPUT='   + str(input.infile),
        'OUTPUT='  + str(output.outfile)
        ])
    shell(command)

# ----------------------------------------------------------------------------- #
rule map_star:
    input:
        infile = rules.sam_to_fastq.outfile,
        genome = rules.make_star_reference.output
    output:
        outfile   = temp(os.path.join(PATH_MAPPED, "{name}", "{genome}","star.Aligned.out.sam"))
    params:
        star       = SOFTWARE['star']['executable'],
        outpth     = os.path.join(PATH_MAPPED, "{name}", "{genome}"),
        threads    = 4,
        mem        = '32G',
        tempdir    = TEMPDIR,
    log:
       log = os.path.join(PATH_LOG, "{name}.{genome}.star.log")

    shell"""
    {param.star} --genomeDir {input.genome} --runThreadN {params.threads} --outFileNamePrefix {params.outpath}/star. --readFilesIn {input.infile}
    """

# ----------------------------------------------------------------------------- #
rule sort_aligned:
    input:
        rules.map_star.outfile
    output:
        outfile   = temp(os.path.join(PATH_MAPPED, "{name}", "{genome}","aligned.sorted.bam"))
    params:
        picard     = SOFTWARE['picard']['executable'],
        java       = SOFTWARE['java']['executable'],
        app_name   = 'SortSam',
        threads    = 1,
        mem        = '5G',
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
        mapped    = rules.sort_aligned.outfile,
        unmapped  = rules.trim_polya.outfile,
        reference = os.path.join(PATH_ANNOTATION, '{genome}', '{genome}.fasta'),
        dict      = rules.fasta_dice.output
    output:
        outfile   = temp(os.path.join(PATH_MAPPED, "{name}", "{genome}","merged.bam"))
    params:
        picard     = SOFTWARE['picard']['executable'],
        java       = SOFTWARE['java']['executable'],
        app_name   = 'MergeBamAlignment',
        threads    = 1,
        mem        = '5G',
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
        infile    = rules.merge_bam.outfile,
        refflat   = rules.gtf_to_refflat.output
    output:
        outfile   = os.path.join(PATH_MAPPED, "{name}", "{genome}","star_gene_exon_tagged.bam")
    params:
        droptools  = SOFTWARE['droptools']['executable'],
        java       = SOFTWARE['java']['executable'],
        app_name   = 'TagReadWithGeneExon',
        threads    = 1,
        mem        = '5G',
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
        'INPUT=f' + str(input.infile),
        'OUTPUT=' + str(output.outfile)
        ])
    shell(command)



# ----------------------------------------------------------------------------- #
rule extract_read_statistics:
    input:
        bamfile = rules.map_scRNA.output.outfile
    output:
        outfile = os.path.join(PATH_MAPPED, "{name}", "{genome}",'{name}_{genome}_ReadStatistics.txt')
    params:
        outname  = "{name}_{genome}",
        threads  = 1,
        mem      = '8G',
        Rscript  = PATH_RSCRIPT
    message: """
            extract_read_statistics:
                input:  {input.bamfile}
                output: {output.outfile}
        """
    run:
        RunRscript(input, output, params, 'Extract_Read_Statistics.R')

# ----------------------------------------------------------------------------- #
rule extract_downstream_statistics:
    input:
        umi_matrix   = rules.get_umi_matrix.output.outfile,
        reads_matrix = rules.get_reads_matrix.output.outfile,
        reads_stats  = rules.extract_read_statistics.output.outfile,
        bam_tag_hist = rules.bam_tag_histogram.output.outfile
    output:
        outfile = os.path.join(PATH_MAPPED, "{name}", "{genome}",'{name}_{genome}_DownstreamStatistics.txt')
    params:
        file_location = os.path.join(PATH_MAPPED, "{name}", "{genome}"),
        outname       = "{name}_{genome}",
        threads       = 1,
        mem           = '8G',
        scriptdir     = SCRIPTDIR,
        Rscript       = RSCRIPT_PATH
    message: """
            extract_downstream_statistics:
                umi:    {input.umi_matrix}
                reads:  {input.reads_matrix}
                stats:  {input.reads_stats}
                output: {output.outfile}
        """
    run:
        RunRscript(input, output, params, 'Extract_Downstream_Statistics.R')


# ----------------------------------------------------------------------------- #
# calculates the number of reads per cell
rule bam_tag_histogram:
    input:
        infile = rules.map_scRNA.output.outfile
    output:
        outfile = os.path.join(PATH_MAPPED, "{name}", "{genome}",'{name}_{genome}_BAMTagHistogram.txt')
    params:
        outdir    = os.path.join(PATH_MAPPED, "{name}", "{genome}"),
        outname   = "{name}_{genome}",
        threads   = 1,
        mem       = '64G'
        java      = SOFTWARE['java']['executable'] ,
        droptools = SOFTWARE['droptools']['executable'],
        tempdir   = TEMPDIR
    message: """
            BamTagHistogram:
                input:  {input.infile}
                output: {output.outfile}
        """
    shell:"""
        {params.java} -XX:ParallelGCThreads={params.threads} -Xmx${params.mem} -Djava.io.tmpdir={params.tempdir} {params.droptools} BAMTagHistogram O={output.outfile} I={input.infile} TAG='XC'
	"""


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
        threads  = 1,
        mem      = '8G',
        cutoff   = 50000,
        Rscript  = RSCRIPT_PATH
    message: """
            find_absolute_read_cutoff:
                input:  {input.infile}
                output: {output.outfile}
        """
    run:
        RunRscript(input, output, params, 'Find_Absolute_Read_Cutoff.R')


# ----------------------------------------------------------------------------- #
# calculates the UMI matrix
rule get_umi_matrix:
    input:
        infile        = rules.map_scRNA.output,
        reads_cutoff  = rules.find_absolute_read_cutoff.output.outfile
    output:
        outfile = os.path.join(PATH_MAPPED, "{name}", "{genome}",'{name}_{genome}_UMI.Matrix.txt')
    params:
        outdir            = os.path.join(PATH_MAPPED, "{name}", "{genome}"),
        outname           = "{name}_{genome}",
        threads           = 1,
        mem               = '8G',
        java              = SOFTWARE['java']['executable'] ,
        droptools         = SOFTWARE['droptools']['executable'],
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


        tool = params.java +
        ' -XX:ParallelGCThreads=' + params.threads +
        ' -Xmx$'                  + params.mem +
        ' -Djava.io.tmpdir='      + params.tempdir +
        params.droptools + ' DigitalExpression'

        command = ' '.join([
        tool,
        'O=' + output.outfile,
        'I=' + str(input.infile),
        'SUMMARY=' + os.path.join(params.outdir, params.outname + '_Summary.txt'),
#         'MIN_NUM_GENES_PER_CELL=' + str(params.genes_per_cell),
        'MIN_NUM_READS_PER_CELL=' + str(reads_cutoff)
        ])
        shell(command)


# ----------------------------------------------------------------------------- #
# calculates the reads matrix - used for PCR duplication estimations
rule get_reads_matrix:
    input:
        infile        = rules.map_scRNA.output,
        reads_cutoff  = rules.find_absolute_read_cutoff.output.outfile
    output:
        outfile = os.path.join(PATH_MAPPED, "{name}", "{genome}",'{name}_{genome}_READS.Matrix.txt')
    params:
        outdir            = os.path.join(PATH_MAPPED, "{name}", "{genome}"),
        outname           = "{name}_{genome}",
        droptools         = SOFTWARE['droptools']['bin'],
        threads           = 1,
        mem               = '8G',
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

        tool = os.path.join(params.droptools,'DigitalExpression')
        command = ' '.join([
        tool,
        'O=' + output.outfile,
        'I=' + str(input.infile),
        'SUMMARY=' + os.path.join(params.outdir, params.outname + '_Summary.txt'),
        'MIN_NUM_READS_PER_CELL=' + str(reads_cutoff),
        'OUTPUT_READS_INSTEAD=T'
        ])
        shell(command)

# ----------------------------------------------------------------------------- #
# convert UMI matrix from txt format into one loom format
rule convert_matrix_from_txt_to_loom:
    input:
        infile        = os.path.join(PATH_MAPPED, "{name}", "{genome}",'{name}_{genome}_UMI.Matrix.txt'),
        gtf           = lambda wildcards: os.path.join(PATH_ANNOTATION, wildcards.genome, '.'.join([wildcards.genome, 'gtf']))
    output:
        outfile       = os.path.join(PATH_MAPPED, "{name}", "{genome}",'{name}_{genome}_UMI.Matrix.loom')
    log:
        log = os.path.join(PATH_LOG, "{name}.{genome}.convert2loom.log")
    message: """
            Convert UMI Matrix from .txt to .loom:
                input:  {input.infile}
                output: {output.outfile}
        """
    shell: "python {PATH_SCRIPT}/convert_matrix_to_loom.py {wildcards.name} {input.infile} {input.gtf} {output.outfile} &> {log.log}"


# ----------------------------------------------------------------------------- #
## combines multiple loom files into one loom file
rule combine_UMI_matrices_into_loom:
    input:
        infile        = lambda wildcards: expand(os.path.join(PATH_MAPPED, "{name}", wildcards.genome, '_'.join(["{name}", wildcards.genome, 'UMI.Matrix.loom'])), name = SAMPLE_NAMES)
    output:
        outfile       = os.path.join(PATH_MAPPED, "{genome}_UMI.loom")
    params:
        tmpfile       = os.path.join(PATH_MAPPED, "{genome}_UMI.loom.tmp")
    log:
        log = os.path.join(PATH_LOG, "{genome}.combine_looms.log")
    message: """
            Combine multiple loom files into one:
                input:  {input.infile}
                output: {output.outfile}
        """
    shell: "echo {input.infile} > {params.tmpfile}; python {PATH_SCRIPT}/combine_UMI_matrices.py {params.tmpfile} {output.outfile}; rm {params.tmpfile} &> {log.log}"


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
        rscript = SOFTWARE['R']['Rscript']
    message: """
            Import loom file, preprocess and save as singleCellExperiment.RDS:
                input:  {input.infile}
                output: {output.outfile}
        """
    shell: "{params.rscript} {PATH_SCRIPT}/convert_loom_to_singleCellExperiment.R --loomFile={input.infile} --metaDataFile={PATH_META_DATA} --genomeBuild={wildcards.genome} --outFile={output.outfile} &> {log.log}"


# ----------------------------------------------------------------------------- #
## Using the preprocessed SingleCellExperiment.RDS file, generates a self-contained HTML report
rule report:
    input:
        infile        = os.path.join(PATH_MAPPED, "{genome}.SingleCellExperiment.RDS")
    output:
        outfile       = os.path.join(PATH_MAPPED, "{genome}.scRNA-Seq.report.html")
    params:
        reportRmd     = os.path.join(PATH_SCRIPT, "scrnaReport.Rmd"),
        rscript       = SOFTWARE['R']['Rscript']
    log:
        log = os.path.join(PATH_LOG, "{genome}.scRNA-Seq.report.log")
    message: """
            Generate an HTML report from SingleCellExperiment.RDS:
                input:  {input.infile}
                output: {output.outfile}
        """
    shell: "{params.rscript} {PATH_SCRIPT}/renderReport.R --reportFile={params.reportRmd} --sceRdsFile={input.infile} --covariates='{COVARIATES}' --prefix={wildcards.genome} --workdir={PATH_MAPPED} &> {log.log}"

# ----------------------------------------------------------------------------- #
rule bam_to_BigWig:
    input:
        bamfile            = rules.map_scRNA.output,
        cell_cutoff_file   = rules.find_absolute_read_cutoff.output.outfile,
        reads_by_cell_file = rules.bam_tag_histogram.output.outfile
    output:
        bwfile = os.path.join(PATH_MAPPED, "{name}", "{genome}",'{name}_{genome}.bw')
    params:
        threads = 1,
        mem     = '16G',
        Rscript = PATH_RSCRIPT
    message: """
            bam_to_BigWig:
                input:  {input.bamfile}
                output: {output.bwfile}
            """
    run:
        RunRscript(input, output, params, 'BamToBigWig.R')


# ----------------------------------------------------------------------------- #
rule fastqc:
    input:
        unpack(get_fastq_files)
    output:
        outfile = os.path.join(PATH_MAPPED, "{name}", "{name}.fastqc.done")
    params:
        outpath = os.path.join(PATH_MAPPED, "{name}"),
        threads = 1,
        mem     = '16G',
        java    = SOFTWARE['java']['executable']
    log:
        log = os.path.join(PATH_LOG, "{name}.fastqc.log")
    message: """
            fastqc:
                input_R1: {input.barcode}
                input_R2: {input.reads}
                output: {output.outfile}
            """
    shell:"""
        fastqc -j {params.java} -t {params.threads} -o {params.outpath} {input.barcode} {input.reads} 2> {log.log}
        touch {output.outfile}
    """
