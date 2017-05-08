

"""
# ---------------------------------------------------------------------------------- #
ssh max-login2
qrsh
cpath='/home/vfranke/bin/Software/miniconda3/bin'
source $cpath/activate p35

SNAKEFILE='/home/vfranke/Projects/FDamm_scRNAseq/Scripts/Snakemake/Snake_Dropseq.py'
WORKDIR='/data/akalin/Projects/FDamm_scRNAseq'
snakemake -R --snakefile $SNAKEFILE --directory $WORKDIR --timestamp --jobs 24 --cluster "qsub -V -l h_vmem={params.mem} -pe smp {params.threads} -l h_rt=36:00:00" --rerun-incomplete --latency-wait 30

# rm $(snakemake --summary | tail -n+2 | cut -f1)

"""


"""
Variables to extract in the config file:
DROPTOOLS PATH
PICARD PATH
STAR PATH
LOCATION OF THE HUMAN AND MOUSE GENOMES
TARGET AND MIXED GENOMES
PROJECT PATH
MARKDOWN PATH


ADDITIONAL WORK:

* scRNA_report - dependency on the get_umi_matrix must be set
* Create rules for making all of the annotation
* The rmarkdown document needs to be updated for the experiment selection
* sort out the getfiles function to be unambiguous
* write checks for data structures
"""


import glob
import os
import re
import subprocess

PICARD = "~/bin/Software/Dropseq/Drop-seq_tools-1.12/3rdParty/picard/picard.jar"
DROPTOOLS="~/bin/Software/Dropseq/Drop-seq_tools-1.12"
STAR="~/bin/Software/Mapping/STAR-STAR_2.4.2a/bin/Linux_x86_64_static/STAR"

BASE = os.path.join(os.environ['BASEPATH'],'Dropseq')
GENOMES = { 'hg19' : BASE + '/hg19/STAR',
# 			'mm9':   BASE + '/mm9/STAR',
            'mm9_hg19': BASE + '/mm9_hg19/STAR'}

REFERENCE = { 'hg19': BASE + '/hg19/hg19.fasta',
# 			  'mm9':  BASE + '/mm9/mm9.fasta',
              'mm9_hg19': BASE + '/mm9_hg19/mm9_hg19.fasta'}

# DEPENDS ON THE RELATIVE PATH - HAVE TO CHANGE THIS FOR BETTER FUNCTIONALITY
# IF I REMOVE THIS THE MARKDOWN DOCUMENT DOES NOT WORK
PROJECT_PATH=os.path.join(os.environ['PROJECTPATH'],'FDamm_scRNAseq')


REPORT_MARKDOWN = '/home/vfranke/Projects/FDamm_scRNAseq/Scripts/Analysis/scRNA_Report.rmd'


BATCH = 'DROP-SEQ-T4-HUMANMOUSE3'


# ----------------------------------------------------------------------------- #
files = [file for file in glob.glob(os.path.join(PROJECT_PATH,"Data",BATCH) + "/*/*.fastq.gz", recursive=True)]
NAMES = list(set([re.sub('_R\\d_.+','',os.path.basename(file)) for file in files]))
MAPPED_DIR = os.path.join(PROJECT_PATH, "Mapped/STAR", BATCH)
QC_PATH = os.path.join(MAPPED_DIR,"QC")

localrules: report_links, mk_outdir

# ----------------------------------------------------------------------------- #
# functions
def getfiles(wc):
	return(dict(zip(['R1','R2'], sorted(list(set(glob.glob('**/' + BATCH + '/**/*' + wc.name + '*.fastq.gz', recursive=True)))))))


# ----------------------------------------------------------------------------- #
REFFLAT = expand(os.path.join(BASE,'{genome}','{genome}.refFlat'), genome=REFERENCE.keys())
DICT    = expand(os.path.join(BASE,'{genome}','{genome}.dict'), genome=REFERENCE.keys())
FASTQC  = expand(os.path.join(MAPPED_DIR, "{name}", 'FastQC', '{name}_R2_001_fastqc.zip'), name=NAMES)
UMI     = expand(os.path.join(MAPPED_DIR, "{name}", "{genome}",'{name}_{genome}_UMI.Matrix.txt'),
name=NAMES, genome=GENOMES.keys())
READ    = expand(os.path.join(MAPPED_DIR, "{name}", "{genome}",'{name}_{genome}_Read.Matrix.txt'), name=NAMES, genome=GENOMES.keys())

REPORT  = expand(os.path.join(MAPPED_DIR, "{name}", "Report",  BATCH + "_" + "{name}.scRNA_Report.html"), name=NAMES)

LINKS   = expand(os.path.join(PROJECT_PATH, "Reports",  BATCH + "_" + "{name}.scRNA_Report.html"), name=NAMES)


# ----------------------------------------------------------------------------- #
rule all:
	input:
		REFFLAT + DICT + FASTQC + UMI + READ + REPORT + LINKS



# ----------------------------------------------------------------------------- #
rule fasta_dict:
	input:
		os.path.join(BASE, '{genome}', '{genome}.fasta')
	output:
		os.path.join(BASE, '{genome}', '{genome}.dict')
	params:
		picard=PICARD,
		threads=1,
		mem='4G'
	message:"""
		Fasta dict...
	"""
	shell:"""
		java -Xmx1200m -jar {params.picard} CreateSequenceDictionary R={input} O={output}
	"""

# ----------------------------------------------------------------------------- #
rule gtf_to_refflat:
	input:
		dict = os.path.join(BASE, '{genome}', '{genome}.dict'),
		gtf = os.path.join(BASE, '{genome}', '{genome}.gtf')
	output:
		os.path.join(BASE, '{genome}', '{genome}.refFlat')
	params:
		droptools=DROPTOOLS,
		threads=1,
		mem='16G'
	message:"""
        refflat...
	"""
	shell:"""
		{params.droptools}/ConvertToRefFlat O={output} ANNOTATIONS_FILE={input.gtf} SEQUENCE_DICTIONARY={input.dict}
	"""

# ----------------------------------------------------------------------------- #
rule fastqc:
	input:
		unpack(getfiles)
	output:
		os.path.join(MAPPED_DIR, "{name}", 'FastQC', '{name}_R2_001_fastqc.zip')
	params:
		outdir = os.path.join(MAPPED_DIR, "{name}", 'FastQC'),
		threads=1,
		mem='4G'
	message:"""
        fastqc...
	"""
	shell:"""
		mkdir {params.outdir}
		fastqc -t 1 -o {params.outdir} {input.R2}
	"""


# ----------------------------------------------------------------------------- #



rule fastq_to_bam:
	input:
		unpack(getfiles)
	output:
		os.path.join(MAPPED_DIR, "{name}", "{name}.fastq.bam")
	params:
		name='{name}',
		picard=PICARD,
		threads=1,
		mem='16G'
	message: """
		Merging Fastq ...
	"""
	shell: """
		java -Xmx12000m  -jar {params.picard} FastqToSam O={output} F1={input.R1} F2={input.R2} QUALITY_FORMAT=Standard SAMPLE_NAME={params.name} SORT_ORDER=queryname
	"""



# ----------------------------------------------------------------------------- #
rule map_dropseq:
	input:
		refflat = rules.gtf_to_refflat.output,
		dict = rules.fasta_dict.output,
		infile = rules.fastq_to_bam.output,
		genome = lambda wc: GENOMES[wc.genome],
		reference = lambda wc: REFERENCE[wc.genome]
	output:
		outfile = os.path.join(MAPPED_DIR, "{name}", "{genome}", "star_gene_exon_tagged.bam")
	params:
		outdir = os.path.join(MAPPED_DIR, "{name}", "{genome}"),
		droptools=DROPTOOLS,
		star = STAR,
		threads=8,
		mem='8G'
	log:
		os.path.join(MAPPED_DIR, "{name}", "{genome}","{name}.log")
	message: """
		Mapping dropseq...
	"""
	shell:"""
		{params.droptools}/Drop-seq_alignment.sh -g {input.genome} -d {params.droptools} -o {params.outdir} -s {params.star} -r {input.reference} {input.infile}
		touch {output.outfile}
	"""


# ----------------------------------------------------------------------------- #
rule get_umi_matrix:
	input:
		infile = rules.map_dropseq.output
	output:
		os.path.join(MAPPED_DIR, "{name}", "{genome}",'{name}_{genome}_UMI.Matrix.txt')
	params:
		outdir = os.path.join(MAPPED_DIR, "{name}", "{genome}"),
		outname = "{name}_{genome}",
		droptools=DROPTOOLS,
		threads=1,
		mem='8G'
	message:"""
		Count UMI ...
	"""
	shell:"""
		{params.droptools}/DigitalExpression O={output} I={input} SUMMARY={params.outdir}/{params.outname}_Summary.txt MIN_NUM_GENES_PER_CELL=10 NUM_CORE_BARCODES=10
	"""

# ----------------------------------------------------------------------------- #
rule get_read_matrix:
	input:
		infile = rules.map_dropseq.output
	output:
		os.path.join(MAPPED_DIR, "{name}", "{genome}",'{name}_{genome}_Read.Matrix.txt')
	params:
		outdir = os.path.join(MAPPED_DIR, "{name}", "{genome}"),
		outname = "{name}_{genome}",
		droptools=DROPTOOLS,
		threads=1,
		mem='8G'
	message:"""
		Count reads ...
	"""
	shell:"""
		{params.droptools}/DigitalExpression O={output} I={input} SUMMARY={params.outdir}/{params.outname}_Summary.txt MIN_NUM_GENES_PER_CELL=10 NUM_CORE_BARCODES=10 OUTPUT_READS_INSTEAD=true
	"""


# ----------------------------------------------------------------------------- #
rule scRNA_report:
	input:
		inpath  =  os.path.join(MAPPED_DIR, "{name}"),
		infiles = [os.path.join(MAPPED_DIR, "{name}", genome,'{name}_'+ genome +'_UMI.Matrix.txt') for genome in GENOMES.keys()]
	output:
		outfile = os.path.join(MAPPED_DIR, "{name}", "Report",  BATCH + "_" + "{name}.scRNA_Report.html")
	params:
		dropseq = BASE,
		report  = REPORT_MARKDOWN,
		threads=1,
		mem='80G'
	message:"""
		scRNA_Report ...
	"""
	shell:
		"""
		Rscript -e 'rmarkdown::render(input = "{params.report}", output_file = "{output.outfile}", knit_root_dir="{input.inpath}", params=list(dropseq="{params.dropseq}"))'
 		"""


# ----------------------------------------------------------------------------- #
rule report_links:
    input:
        infile = rules.scRNA_report.output
    output:
        os.path.join(PROJECT_PATH, "Reports", BATCH + "_" + "{name}.scRNA_Report.html")
    params:
        threads=1,
        mem='1G'
    message:"""
        Report links ...
	"""
	shell:"""
    	ln -s {input.infile} {output}
	"""
