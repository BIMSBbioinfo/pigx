"""
jump conda3
source ./activate p35
SNAKEFILE='/home/vfranke/Projects/AAkalin_PIX/Scripts/ChIP/Snake_ChIPseq.py'
WORKDIR='/data/local/vfranke/AAkalin_PIX/ChIP'
CONFIGFILE='/home/vfranke/Projects/AAkalin_PIX/Scripts/ChIP/config.yaml'

# Beast run
snakemake -R --snakefile $SNAKEFILE --directory $WORKDIR --jobs 4 --rerun-incomplete --configfile $CONFIGFILE --latency-wait 30

snakemake -R --snakefile $SNAKEFILE --directory $WORKDIR --jobs 4 --rerun-incomplete --configfile $CONFIGFILE --latency-wait 30 --dryrun

# Visualize the graph horizontally
snakemake --dag --snakefile $SNAKEFILE --directory $WORKDIR --configfile $CONFIGFILE | perl -pe 's|graph\[|graph[rankdir=LR,|'| dot -T png > Snake_ChIPseq.png

# Cluster run
snakemake -R --snakefile $SNAKEFILE --directory $WORKDIR --timestamp --jobs 24 --cluster "qsub -V -l h_vmem={params.mem} -pe smp {params.threads} -l h_rt=36:00:00" --rerun-incomplete --latency-wait 30


rm $(snakemake --summary | tail -n+2 | cut -f1)

# TODO
* 0. Make test data
* 1. Variables into config file:
*    - genome path
*    - read extension

* 2. Automatic reference generation
* 3. check globbing
* 4. Add peak calling
5. Add peak QC
6. Add interactive macs parameters
7. Add tests for rules
8. Add checks for the config file
9. Add possible pseudonims for samples
10. Extent to paired end reads
    - Testing for paired end reads
    - Mapping
    - Fastqc for paired end reads

11. Run multiqc
12. Add samplesheet To yaml file
13. Check for params:idr
* 14. Extract the input file for the fastqs as a config parameter
* 15. application placeholders in config files
16. Format messages

17. cofigure the pipeline for broad histone data

Integrate_Expression_Mouse_Hamster_Difference
1. Figure out how to prevent bombing

### 13. config file check

Q:
1. How to write dynamic rules? - that the rules execute based on the input parameters?
    - if the user supplies the genome index, do not make the index
2. What to do downstream? Some basic plots? ChromHMM
3. How to speed up the extension?

"""

SCRIPT_PATH = 'Scripts'
PARAMS_PATH = '.'

# ---------------------------------------------------------------------------- #
import glob
import os
import re
import sys
import yaml

from SnakeFunctions import *
localrules: makelinks


# ---------------------------------------------------------------------------- #
# check config validity
if check_config(config) == 1:
    print('config.yaml is not properly formatted - exiting')
    quit()

# ---------------------------------------------------------------------------- #
# Variable definition
GENOME       = config['genome']['name']
GENOME_FASTA = config['genome']['fasta']
NAMES        = config['samples'].keys()
PATH_FASTQ   = config['fastq']
SOFTWARE     = config['software']
PARAMS 		 = config['params']

# Directory structure definition
PATH_MAPPED = "Mapped/Bowtie"
PATH_QC     = "FastQC"
PATH_INDEX  = 'Bowtie_Index'
PATH_LOG    = 'Log'
PATH_PEAK   = 'Peaks/MACS2'
PATH_BW     = 'Links_bw'
PATH_IDR    = 'IDR'

# ---------------------------------------------------------------------------- #
# Function parameter
with open(os.path.join(PARAMS_PATH,'app_params.yaml'), 'r') as stream:
    app_params = yaml.load(stream)


# ---------------------------------------------------------------------------- #
# config defaults
if not ('extend' in config['params'].keys()):
    config['params']['extend'] = 0


if GENOME_FASTA == None:
	prefix_default = ''
	GENOME_FASTA = ''
else:
	prefix_default = os.path.join(PATH_INDEX, GENOME_NAME)


PREFIX = os.path.join(set_default('index', prefix_default, config), GENOME)



# ----------------------------------------------------------------------------- #
# RULE ALL
INDEX   = [PREFIX + '.1.ebwt']
MAPPING = expand(os.path.join(PATH_MAPPED, "{name}", "{name}.sorted.bam.bai"), name=NAMES)
FASTQC  = expand(os.path.join(PATH_QC,     "{name}", "{name}.fastqc.done"), name=NAMES)
BW      = expand(os.path.join(os.getcwd(), PATH_MAPPED, "{name}", "{name}.bw"),  name=NAMES)
LINKS   = expand(os.path.join(PATH_BW,  "{ex_name}.bw"),  ex_name=NAMES)
MACS    = expand(os.path.join(PATH_PEAK,    "{name}", "{name}_peaks.narrowPeak"), name=config['peak_calling'].keys())
IDR     = expand(os.path.join(PATH_IDR,    "{name}", "{name}.narrowPeak"), name=config['idr'].keys())

rule all:
    input:
        INDEX + MAPPING + BW + LINKS + FASTQC + MACS + IDR

# ----------------------------------------------------------------------------- #
rule bowtie_build:
    input:
        genome = GENOME_FASTA
    output:
        outfile = PREFIX + '.1.ebwt'
    params:
        prefix = PREFIX,
        threads = 1,
        mem = '32G',
        bowtie_bild = SOFTWARE['bowtie-build']
    log:
        os.path.join(PATH_LOG, "bowtie_build.log")
    message:
    	"""
        	Constructing bowtie index:
        		input : {input.genome}
        		output: {output.outfile}
    	"""
    shell:"""
        {params.bowtie_build} {input.genome} {params.prefix} 2> {log}
    """


#----------------------------------------------------------------------------- #
def fastq_input(wc):
    samps = config['samples'][wc.name]

    if type(samps) is str:
        samps = [samps]

    infiles = [os.path.join(PATH_FASTQ, i) for i in samps]
    return(infiles)

rule bowtie:
    input:
        infile = fastq_input,
        genome = rules.bowtie_build.output.outfile
    output:
        bamfile = os.path.join(PATH_MAPPED, "{name}/{name}.bam")
    params:
        threads=8,
        mem='16G',
        bowtie = SOFTWARE['bowtie'],
        samtools = SOFTWARE['samtools']
    log:
        log = os.path.join(PATH_MAPPED, "{name}/{name}.bowtie.log")
    message:
    """
        Mapping with bowtie:
        	sample: {input.infile}
        	genome: {input.genome}
        	output: {output.bamfile}
    """
    run:
        genome = input.genome.replace('.1.ebwt','')

        command = input.infile[0] + ' | ' + params.bowtie + ' -p ' + str(params.threads) +' -S -k 1 -m 1 --best --strata ' + genome  + ' -  2> ' + log.log + ' | '+ params.samtools +' view -bhS > ' + output.bamfile

        if os.path.splitext(input.infile[0])[1] == 'gz':
            command = 'zcat ' + command
        else:
            command = 'cat ' + command

        shell(command)


# ----------------------------------------------------------------------------- #
rule samtools_sort:
    input:
        os.path.join(PATH_MAPPED, "{name}/{name}.bam")
    output:
        os.path.join(PATH_MAPPED, "{name}/{name}.sorted.bam")
    params:
        threads = 4,
        mem = '16G',
        samtools = SOFTWARE['samtools']
    message:
    	"""
        	Sorting mapped reads:
        		input: {input}
        		output: {output}
	    """
    shell: """
        {params.samtools} sort --threads {params.threads} -o {output} {input}
    """

# ----------------------------------------------------------------------------- #
rule samtools_index:
    input:
        os.path.join(PATH_MAPPED, "{name}/{name}.sorted.bam")
    output:
        os.path.join(PATH_MAPPED, "{name}/{name}.sorted.bam.bai")
    params:
        threads = 1,
        mem = '8G',
        samtools = SOFTWARE['samtools']
    message:
    """
        Indexing bam file:\n
        	input: {input}
    """
    shell:
        "{params.samtools} index {input}"


# ----------------------------------------------------------------------------- #
rule fastqc:
    input:
        infile = fastq_input
    output:
        outfile = os.path.join(PATH_QC, "{name}", "{name}.fastqc.done")
    params:
        outpath = os.path.join(PATH_QC, "{name}"),
        threads = 1,
        mem = '8G',
        fastqc = SOFTWARE['fastqc']
    message:
        """
        	FastQC:
        		input: {input.infile}
        		output: {output.outfile}
        """
    shell: """
        {params.fastqc} --outdir {params.outpath} --extract -f fastq -t 1 {input.infile}
        touch {output.outfile}
    """

# ----------------------------------------------------------------------------- #
rule chrlen:
    input:
        os.path.join(PATH_MAPPED, "{name}/{name}.sorted.bam")
    output:
        os.path.join(PATH_MAPPED, "{name}/chrlen.txt")
    params:
        threads = 1,
        mem = '4G',
        samtools = SOFTWARE['samtools'],
        perl = SOFTWARE['perl']
    shell:"""
        {params.samtools} view -H {input} | grep @SQ | {params.perl} -pe 's/^@.+?://;s/LN://' > {output}
    """

rule bam2bed:
    input:
        file = os.path.join(PATH_MAPPED, "{name}/{name}.sorted.bam"),
        chrlen = rules.chrlen.output
    output:
        os.path.join(PATH_MAPPED, "{name}/{name}.bed")
    params:
        extend=config['params']['extend'],
        threads = 1,
        mem = '16G',
        bamToBed = SOFTWARE['bamToBed']
    shell: """
        {params.bamToBed} -i {input.file} > {output}
    """

rule bam2bed_extend:
    input:
        file = rules.bam2bed.output,
        chrlen = rules.chrlen.output
    output:
        file = os.path.join(PATH_MAPPED, "{name}/{name}.ext.bed")
    params:
        extend = config['params']['extend'],
        threads = 1,
        mem = '16G',
    message:
    	"""
    		Extending reads for coverage file:
    			input : {input.file}
    			extend: {params.extend}
    			output: {output.file}
	    """
    script:
    	os.path.join(SCRIPT_PATH, 'Extend_Regions.R')

#
rule bam2bigWig:
    input:
        chrlen = rules.chrlen.output,
        bedfile = rules.bam2bed_extend.output
    output:
        os.path.join(os.getcwd(), PATH_MAPPED, "{name}/{name}.bw")
    params:
        threads = 1,
        mem = '16G',
        bedgraph = os.path.join(PATH_MAPPED, "{name}/{name}.bedGraph"),
        genomeCoverageBed = SOFTWARE['genomeCoverageBed'],
        wigToBigWig = SOFTWARE['wigToBigWig']
    message:
    """
    	Making bigWig:
    		input : {input.bedfile}
    		output: {output}
    """
    shell: """
        {params.genomeCoverageBed} -i {input.bedfile} -g {input.chrlen} -bg > {params.bedgraph}
        {params.wigToBigWig} {params.bedgraph} {input.chrlen} {output}
    """

# ----------------------------------------------------------------------------- #
rule makelinks:
    input:
        file = os.path.join(os.getcwd(), PATH_MAPPED, "{name}", "{name}" + '.bw')
    output:
        os.path.join(PATH_BW, "{name}.bw")
    shell: """
        ln -s {input.file} {output}
    """


# ----------------------------------------------------------------------------- #
def get_sample_macs(wc):
    name = config['peak_calling'][wc.name]
    samps = dict(zip(name.keys(),[os.path.join('Mapped','Bowtie',i, i + '.sorted.bam') for i in name.values()]))
    return(samps)

rule macs2:
    input:
        unpack(get_sample_macs)
    output:
    	outfile = os.path.join(PATH_PEAK, "{name}", "{name}_peaks.narrowPeak")
    params:
    	outpath = os.path.join(PATH_PEAK, "{name}"),
        name = "{name}",
        threads = 1,
        mem = '16G',
        macs2 = SOFTWARE['macs2']
        params = PARAMS['macs2']
    log:
    	os.path.join(PATH_LOG, 'macs.log')
    message:
    """
    	Running macs2:
			ChIP:   {input.ChIP}
			Cont:   {input.Cont}
			output: {output.outfile}
    """
    run:
        command = " ".join(
        [params.macs2, 'callpeak',
        '-t', input.ChIP,
        '-c', input.Cont,
        '--outdir', params.outpath,
        '-n', params.name,
        join_params(params.macs2, app_params, config),
        '2>', log
        ])
    	shell(command)


# ----------------------------------------------------------------------------- #
def get_sample_idr(wc):
    name = config['idr'][wc.name]
    samps = dict(zip(name.keys(),[os.path.join(PATH_PEAK, i, i + '_peaks.narrowPeak') for i in name.values()]))
    return(samps)

rule idr:
	input:
		unpack(get_sample_idr)
	output:
		outfile = os.path.join(PATH_IDR, "{name}", "{name}.narrowPeak")
	params:
		threads = 1,
		mem = '8G',
		idr = SOFTWARE['idr']
	log:
		os.path.join(PATH_LOG, 'idr.log')
	message:
		"""
			Running IDR2:
				input : {input.ChIP1} {input.ChIP2}
				output: {output.outfile}
		"""
	run:
		command = " ".join(
        [params.idr,
        '--samples', input.ChIP1, input.ChIP2,
        '--input-file-type', 'narrowPeak',
        '--rank', 'q.value',
        '--output-file', output.outfile,
        '-l', log,
        '--plot',
        join_params(params.idr, app_params, config)
        ])
        shell(command)
