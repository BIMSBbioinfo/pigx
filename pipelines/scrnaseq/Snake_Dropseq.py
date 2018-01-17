"""
# ---------------------------------------------------------------------------------- #
SNAKEFILE='/home/vfranke/Projects/AAkalin_PIX_scRNA/Snake_Dropseq.py'
WORKDIR='/data/local/vfranke/AAkalin_PIX/scRNA'
CONFIGFILE='/home/vfranke/Projects/AAkalin_PIX_scRNA/Config_scRNA.yaml'
PATH='/home/vfranke/bin/Software/miniconda3/envs/p35/bin:/usr/local/bin:/usr/bin:/bin:/home/vfranke/.guix-profile/bin:/home/vfranke/.guix-profile/sbin:/home/vfranke/bin'
# Beast run

snakemake -R --snakefile $SNAKEFILE --directory $WORKDIR --jobs 12 --rerun-incomplete --configfile $CONFIGFILE --latency-wait 30 --dryrun

snakemake -R --snakefile $SNAKEFILE --directory $WORKDIR --jobs 4 --rerun-incomplete --configfile $CONFIGFILE --latency-wait 30
"""


"""

TODO:


Variables to extract in the config file:
DROPTOOLS
PICARD
STAR
PROJECT PATH
MARKDOWN PATH


"""

# ----------------------------------------------------------------------------- #
# libraries and constants
import glob
import os
import re
import subprocess
import yaml

localrules: report_links, mk_outdir

PATH_SCRIPT = os.path.join(workflow.basedir,'Scripts')
PATH_RULES  = os.path.join(workflow.basedir,'Rules')

# ----------------------------------------------------------------------------- #
# Software parameters
with open(os.path.join(workflow.basedir, 'software_scRNA.yaml'), 'r') as stream:
    SOFTWARE_CONFIG = yaml.load(stream)

# Function parameter
APP_PARAMS = SOFTWARE_CONFIG['params']
SOFTWARE = SOFTWARE_CONFIG['software']
# ----------------------------------------------------------------------------- #
# Variables
GENOME_NAME_PRIMARY = config['annotation']['primary']['genome']['name']
REFERENCE_NAMES = [GENOME_NAME_PRIMARY]

SAMPLE_NAMES = config['samples'].keys()
PARAMS = config['params']

# ----------------------------------------------------------------------------- #
# PATHS
PATH_FASTQ = config['fastq']
PATH_ANNOTATION = 'Annotation'
PATH_ANNOTATION_PRIMARY = os.path.join(PATH_ANNOTATION, GENOME_NAME_PRIMARY)
PATH_REFERENCE_PRIMARY  = config['annotation']['primary']['genome']['fasta']
PATH_GTF_PRIMARY        = config['annotation']['primary']['gtf']


GENOME_NAME_MIX = None
if 'secondary' in set(config['annotation'].keys()):
    GENOME_NAME_SECONDARY    = config['annotation']['secondary']['genome']['name']
    PATH_REFERENCE_SECONDARY = config['annotation']['secondary']['genome']['fasta']
    PATH_GTF_SECONDARY       = config['annotation']['secondary']['gtf']

    GENOME_NAME_MIX = GENOME_NAME_PRIMARY + '_' + GENOME_NAME_SECONDARY
    PATH_ANNOTATION_MIX = os.path.join(PATH_ANNOTATION, GENOME_NAME_MIX)
    REFERENCE_NAMES = REFERENCE_NAMES + [GENOME_NAME_MIX]


PATH_BASE = os.path.join(os.environ['BASEPATH'],'Dropseq')

PATH_MAPPED   = 'Mapped'
PATH_LOG      = 'Log'
PATH_MARKDOWN = os.path.join(workflow.basedir,'Templates/scRNA_Report.rmd')


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

# ----------------------------------------------------------------------------- #
# Bam To BigWig
BIGWIG = expand(os.path.join(PATH_MAPPED, "{name}", "{genome}",'{name}_{genome}.bw'), genome = REFERENCE_NAMES, name = SAMPLE_NAMES)


# ----------------------------------------------------------------------------- #
RULE_ALL = []
RULE_ALL = RULE_ALL + [LINK_REFERENCE_PRIMARY, LINK_GTF_PRIMARY]

if len(COMBINE_REFERENCE) > 0:
    RULE_ALL = RULE_ALL + COMBINE_REFERENCE

RULE_ALL = RULE_ALL + DICT + REFFLAT + MAKE_STAR_INDEX + FASTQC + MERGE_FASTQ_TO_BAM + MAP_scRNA + READ_STATISTICS + BAM_HISTOGRAM + FIND_READ_CUTOFF + UMI + BIGWIG
print(RULE_ALL)


# FASTQC  = expand(os.path.join(MAPPED_DIR, "{name}", 'FastQC', '{name}_R2_001_fastqc.zip'), name=NAMES)
# UMI     = expand(os.path.join(MAPPED_DIR, "{name}", "{genome}",'{name}_{genome}_UMI.Matrix.txt'),
# name=NAMES, genome=GENOMES.keys())
# READ    = expand(os.path.join(MAPPED_DIR, "{name}", "{genome}",'{name}_{genome}_Read.Matrix.txt'), name=NAMES, genome=GENOMES.keys())
#
# REPORT  = expand(os.path.join(MAPPED_DIR, "{name}", "Report",  BATCH + "_" + "{name}.scRNA_Report.html"), name=NAMES)
#
# LINKS   = expand(os.path.join(PROJECT_PATH, "Reports",  BATCH + "_" + "{name}.scRNA_Report.html"), name=NAMES)
# 		REFFLAT + DICT + FASTQC + UMI + READ + REPORT + LINKS


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
        picard=SOFTWARE['picard'],
        threads=1,
        mem='4G'
    log:
        os.path.join(PATH_LOG, '{genome}.fasta_dict.log')
    message:
        """
            Fasta dict:
                input  : {input}
                output : {output}
        """
    shell:"""
        java -Xmx1200m -jar {params.picard} CreateSequenceDictionary R={input} O={output} 2> {log}
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
        threads=1,
        mem='4G'
    message:
        """
            Changing GTF id:
                input  : {input}
                output : {output}
        """
    script:
        os.path.join(PATH_SCRIPT, 'change_gtf_id.R')



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
        droptools=SOFTWARE['droptools'],
        threads=1,
        mem='50G'
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
        {params.droptools}/ConvertToRefFlat O={output} ANNOTATIONS_FILE={input.gtf} SEQUENCE_DICTIONARY={input.dict} 2> {log}
    """

# # ----------------------------------------------------------------------------- #
def get_fastq_files(wc):

    h = {'barcode' : os.path.join(PATH_FASTQ, config['samples'][wc.name]['barcode']),
         'reads'   : os.path.join(PATH_FASTQ, config['samples'][wc.name]['reads'])}
    return h

rule merge_fastq_to_bam:
    input:
        unpack(get_fastq_files)
    output:
        os.path.join(PATH_MAPPED, "{name}", "{name}.fastq.bam")
    params:
        name='{name}',
        picard=SOFTWARE['picard'],
        threads=1,
        mem='16G'

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
    java -Xmx12000m  -jar {params.picard} FastqToSam O={output} F1={input.barcode} F2={input.reads} QUALITY_FORMAT=Standard SAMPLE_NAME={params.name} SORT_ORDER=queryname 2> {log}
    """


# ----------------------------------------------------------------------------- #
# Maps reads using Dropseq tools and STAR
rule map_scRNA:
    input:
        refflat   = rules.gtf_to_refflat.output,
        dict      = rules.fasta_dict.output,
        index     = rules.make_star_reference.output,
        infile    = rules.merge_fastq_to_bam.output,
        reference = os.path.join(PATH_ANNOTATION, '{genome}', '{genome}.fasta'),
    output:
        outfile   = os.path.join(PATH_MAPPED, "{name}", "{genome}","star_gene_exon_tagged.bam")
    params:
        outdir    = os.path.join(PATH_MAPPED, "{name}", "{genome}"),
        genome    = os.path.join(PATH_ANNOTATION, '{genome}', 'STAR_INDEX'),
        droptools = SOFTWARE['droptools'],
        star      = SOFTWARE['star'],
        threads   = 8,
        mem       = '35G'
    log:
       log = os.path.join(PATH_LOG, "{name}.{genome}.STAR.log")
    message: """
            Mapping scRNA:
                input:
                    dict   : {input.dict}
                    reads  : {input.infile}
                    genome : {params.genome}
        """
    shell:"""
        {params.droptools}/Drop-seq_alignment.sh -g {params.genome} -d {params.droptools} -o {params.outdir} -s {params.star} -r {input.reference} {input.infile} 2> {log}
        touch {output.outfile}
    """

# ----------------------------------------------------------------------------- #
rule extract_read_statistics:
    input:
        bamfile = rules.map_scRNA.output.outfile
    output:
        outfile = os.path.join(PATH_MAPPED, "{name}", "{genome}",'{name}_{genome}_ReadStatistics.txt')
    params:
        outname  = "{name}_{genome}",
        threads  = 1,
        mem      = '8G'
    message: """
            extract_read_statistics:
                input:  {input.bamfile}
                output: {output.outfile}
        """
    script:
        os.path.join(PATH_SCRIPT, 'Extract_Read_Statistics.R')

# ----------------------------------------------------------------------------- #
# calculates the number of reads per cell
rule bam_tag_histogram:
    input:
        infile = rules.map_scRNA.output.outfile
    output:
        outfile = os.path.join(PATH_MAPPED, "{name}", "{genome}",'{name}_{genome}_BAMTagHistogram.txt')
    params:
        outdir            = os.path.join(PATH_MAPPED, "{name}", "{genome}"),
        outname           = "{name}_{genome}",
        droptools         = SOFTWARE['droptools'],
        threads           = 1,
        mem               = '64G'
    message: """
            BamTagHistogram:
                input:  {input.infile}
                output: {output.outfile}
        """
    shell:"""
        {params.droptools}/BAMTagHistogram O={output.outfile} I={input.infile} TAG='XC'
		"""

# ----------------------------------------------------------------------------- #

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
        cutoff   = 50000
    message: """
            find_absolute_read_cutoff:
                input:  {input.infile}
                output: {output.outfile}
        """
    script:
        os.path.join(PATH_SCRIPT, 'Find_Absolute_Read_Cutoff.R')
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
        droptools         = SOFTWARE['droptools'],
        threads           = 1,
        mem               = '8G',
        genes_per_cell    = PARAMS['genes_per_cell'],
        num_core_barcodes = PARAMS['num_core_barcodes']
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
#         'MIN_NUM_GENES_PER_CELL=' + str(params.genes_per_cell),
        'MIN_NUM_READS_PER_CELL=' + str(reads_cutoff)
        ])
        shell(command)


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
        mem     = '16G'
    message: """
            bam_to_BigWig:
                input:  {input.bamfile}
                output: {output.bwfile}
            """
    script:
        os.path.join(PATH_SCRIPT, 'BamToBigWig.R')


# ----------------------------------------------------------------------------- #
rule fastqc:
    input:
        unpack(get_fastq_files)
    output:
        outfile = os.path.join(PATH_MAPPED, "{name}","{name}.fastqc.done")
    params:
        outpath = os.path.join(PATH_MAPPED, "{name}"),
        threads = 1,
        mem     = '16G',
        java    = SOFTWARE['java']
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




# ----------------------------------------------------------------------------- #
# rule fastqc:
# 	input:
# 		unpack(getfiles)
# 	output:
# 		os.path.join(MAPPED_DIR, "{name}", 'FastQC', '{name}_R2_001_fastqc.zip')
# 	params:
# 		outdir = os.path.join(MAPPED_DIR, "{name}", 'FastQC'),
# 		threads=1,
# 		mem='4G'
# 	message:"""
#         fastqc...
# 	"""
# 	shell:"""
# 		mkdir {params.outdir}
# 		fastqc -t 1 -o {params.outdir} {input.R2}
# 	"""
#
#

#
#
#

#

# 	"""
#
# # ----------------------------------------------------------------------------- #
# rule get_read_matrix:
# 	input:
# 		infile = rules.map_dropseq.output
# 	output:
# 		os.path.join(MAPPED_DIR, "{name}", "{genome}",'{name}_{genome}_Read.Matrix.txt')
# 	params:
# 		outdir = os.path.join(MAPPED_DIR, "{name}", "{genome}"),
# 		outname = "{name}_{genome}",
# 		droptools=DROPTOOLS,
# 		threads=1,
# 		mem='8G'
# 	message:"""
# 		Count reads ...
# 	"""
# 	shell:"""
# 		{params.droptools}/DigitalExpression O={output} I={input} SUMMARY={params.outdir}/{params.outname}_Summary.txt MIN_NUM_GENES_PER_CELL=10 NUM_CORE_BARCODES=10 OUTPUT_READS_INSTEAD=true
# 	"""
#
#
# # ----------------------------------------------------------------------------- #
# rule scRNA_report:
# 	input:
# 		inpath  =  os.path.join(MAPPED_DIR, "{name}"),
# 		infiles = [os.path.join(MAPPED_DIR, "{name}", genome,'{name}_'+ genome +'_UMI.Matrix.txt') for genome in GENOMES.keys()]
# 	output:
# 		outfile = os.path.join(MAPPED_DIR, "{name}", "Report",  BATCH + "_" + "{name}.scRNA_Report.html")
# 	params:
# 		dropseq = BASE,
# 		report  = REPORT_MARKDOWN,
# 		threads=1,
# 		mem='80G'
# 	message:"""
# 		scRNA_Report ...
# 	"""
# 	shell:
# 		"""
# 		Rscript -e 'rmarkdown::render(input = "{params.report}", output_file = "{output.outfile}", knit_root_dir="{input.inpath}", params=list(dropseq="{params.dropseq}"))'
#  		"""
#
#
# # ----------------------------------------------------------------------------- #
# rule report_links:
#     input:
#         infile = rules.scRNA_report.output
#     output:
#         os.path.join(PROJECT_PATH, "Reports", BATCH + "_" + "{name}.scRNA_Report.html")
#     params:
#         threads=1,
#         mem='1G'
#     message:"""
#         Report links ...
# 	"""
# 	shell:"""
#     	ln -s {input.infile} {output}
# 	"""
