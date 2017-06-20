"""
jump conda3
source ./activate p35
SNAKEFILE='/home/vfranke/Projects/AAkalin_PIX/Snake_ChIPseq.py'
WORKDIR='/data/akalin/vfranke/AAkalin_PIX/ChIP'
CONFIGFILE='/home/vfranke/Projects/AAkalin_PIX/config.yaml'

# Beast run
snakemake -R --snakefile $SNAKEFILE --directory $WORKDIR --jobs 4 --rerun-incomplete --configfile $CONFIGFILE --latency-wait 30

snakemake -R --snakefile $SNAKEFILE --directory $WORKDIR --jobs 4 --rerun-incomplete --configfile $CONFIGFILE --latency-wait 30 --dryrun

# Visualize the graph horizontally
snakemake --dag --snakefile $SNAKEFILE --directory $WORKDIR --configfile $CONFIGFILE | perl -pe 's|graph\[|graph[rankdir=LR,|'| dot -T png > Snake_ChIPseq.png

# Cluster run
snakemake -R --snakefile $SNAKEFILE --directory $WORKDIR --timestamp --jobs 24 --cluster "qsub -V -l h_vmem={params.mem} -pe smp {params.threads} -l h_rt=36:00:00" --rerun-incomplete --latency-wait 30

for i in `ls $WORKDIR | grep -v Fastq$ | grep -v Genome`;do rm -r $WORKDIR/$i;done

rm $(snakemake --snakefile $SNAKEFILE --directory $WORKDIR --configfile $CONFIGFILE --summary | tail -n+2 | cut -f1)

# TODO
## Main
[#]. enable sample specific parameter deposition - works for macs
[#]. write tests
[]. set default app parameters in config
    - works for macs2
[]. the pipeline does not work without app params - change config check    
[]. add check for gzipped bowtie reference
[]. Add samplesheet To yaml file
[]. write markdown
    
[]. write yaml schema
    https://github.com/Grokzen/pykwalify,
    http://www.kuwata-lab.com/kwalify/ruby/users-guide.01.html#schema

[]. extension for paired end data to be automatically determined from the pairs
[]. add scalling to the bedgraph construction
[]. make BigWigExtend a streaming function
    
[]. Tests:
        - for no control sample
        - for multiple chip and control samples
        - add test if control not set in yaml file - check for multiple types of input
        - check whether IDR test works

[]. Tests for hub

## Additional
[] Add peak QC
[]. Run multiqc
[]. cofigure the pipeline for broad histone data
[]. Configure the pipeline for differential peak calling
[]. enable trimming
[]. rewrite fastqc to go through each individual file

[]. Automatic UCSC hub setup
[]. Add test genome data (genome subset)
[]. Set default for genome name if genome specified but genome name is not
[]. Delete temporary files - bed and bedGraph files
[]. Sample specific read extension


DONE
[*] 0. Make test data
[*] 1. Variables into config file:
[*]    - genome path
[*]    - read extension
[*] 2. Automatic reference generation
[*] 3. check globbing
[*] 4. Add peak calling
[*]. Extract the input file for the fastqs as a config parameter
[*]. application placeholders in config files
[*]. Format messages
[*] Add interactive macs parameters
[*] Add possible pseudonims for samples

172405
[*]. write checks for the config proper formatting :
    - ChIP and Cont parameters in peak calling
    - idr samples
[*]. Check for params:idr
[*] Add paired end test data
[*]. Extent to paired end reads
    - Testing for paired end reads
    - Mapping
    - Fastqc for paired end reads
    
170614
[*]. extended the pipeline to accept multiple ChIP/Cont samples for peak calling
[*].  add support for peak calling without control:
    peak calling can be run without control samples - useful for ATAC data

170615    
[*]. Enable running subsections of the pipeline:
    IDR and PEAK calling are not obligatory




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
import fnmatch
import os
import re
import sys
import yaml

from SnakeFunctions import *
from Check_Config import *
localrules: makelinks


# ---------------------------------------------------------------------------- #
# check config validity
if check_config(config) == 1:
    quit()


# ---------------------------------------------------------------------------- #
# Software parameters
with open(os.path.join(workflow.basedir, 'software.yaml'), 'r') as stream:
    SOFTWARE_CONFIG = yaml.load(stream)

# Function parameter
APP_PARAMS = SOFTWARE_CONFIG['params']
# ---------------------------------------------------------------------------- #
# Variable definition
GENOME       = config['genome']['name']
GENOME_FASTA = config['genome']['fasta']
NAMES        = config['samples'].keys()
PATH_FASTQ   = config['fastq']
PARAMS       = config['params']


SOFTWARE     = SOFTWARE_CONFIG['software']


# Directory structure definition
PATH_MAPPED = "Mapped/Bowtie"
PATH_QC     = "FastQC"
PATH_INDEX  = 'Bowtie2_Index'
PATH_LOG    = 'Log'
PATH_PEAK   = 'Peaks/MACS2'
PATH_BW     = 'BigWig'
PATH_IDR    = 'IDR'
PATH_HUB    = 'UCSC_HUB'

# Hub variables which describe the types of files that can be used in the hub
TRACK_PATHS = {
    'bigWig' : {'path': PATH_MAPPED, 'suffix': 'bw', 'type':'bigWig'},
    'macs' :{'path': PATH_PEAK, 'suffix': 'bb', 'type':'bigBed'},
    'idr' : {'path': PATH_IDR, 'suffix': 'IDR.narrowPeak'}
}


# ---------------------------------------------------------------------------- #
# config defaults
if not ('extend' in config['params'].keys()):
    config['params']['extend'] = 0


if GENOME_FASTA == None:
    prefix_default = ''
    GENOME_FASTA = ''
else:
    prefix_default = os.path.join(PATH_INDEX, GENOME)


PREFIX = os.path.join(set_default('index', prefix_default, config), GENOME)
print(PREFIX)


# ----------------------------------------------------------------------------- #
# RULE ALL
COMMAND = []
INDEX   = [PREFIX + '.1.bt2']
BOWTIE2 = expand(os.path.join(PATH_MAPPED, "{name}", "{name}.sorted.bam.bai"), name=NAMES)
FASTQC  = expand(os.path.join(PATH_QC,     "{name}", "{name}.fastqc.done"), name=NAMES)
BW      = expand(os.path.join(os.getcwd(), PATH_MAPPED, "{name}", "{name}.bw"),  name=NAMES)
LINKS   = expand(os.path.join(PATH_BW,  "{ex_name}.bw"),  ex_name=NAMES)

COMMAND = COMMAND + INDEX + BOWTIE2 + BW + LINKS + FASTQC
# ----------------------------------------------------------------------------- #
if 'peak_calling' in set(config.keys()):
    MACS    = expand(os.path.join(PATH_PEAK,  "{name}", "{name}_peaks.narrowPeak"),     name=config['peak_calling'].keys())
    QSORT   = expand(os.path.join(PATH_PEAK,  "{name}", "{name}_qsort.narrowPeak"), name=config['peak_calling'].keys())
    COMMAND = COMMAND + MACS + QSORT 

# ----------------------------------------------------------------------------- #
if 'idr' in set(config.keys()):
    IDR     = expand(os.path.join(PATH_IDR,    "{name}", "{name}.narrowPeak"),     name=config['idr'].keys())
    COMMAND = COMMAND + IDR

# ----------------------------------------------------------------------------- #
HUB_NAME = None
if 'hub' in set(config.keys()):
    HUB_NAME = config['hub']['name']
    HUB = [os.path.join(PATH_HUB, HUB_NAME, 'done.txt')]
    BB  = expand(os.path.join(PATH_PEAK,  "{name}", "{name}.bb"),  name=config['peak_calling'].keys())
    COMMAND = COMMAND + BB + HUB 

rule all:
    input:
        COMMAND

# ----------------------------------------------------------------------------- #
rule bowtie2_build:
    input:
        genome = GENOME_FASTA
    output:
        outfile = PREFIX + '.1.bt2'
    params:
        prefix = PREFIX,
        threads = 1,
        mem = '32G',
        bowtie2_build = SOFTWARE['bowtie2-build']
    log:
        os.path.join(PATH_LOG, "bowtie2_build.log")
    message:
        """
            Constructing bowtie2 index:
                input : {input.genome}
                output: {output.outfile}
        """
    shell:"""
        {params.bowtie2_build} {input.genome} {params.prefix} 2> {log}
    """


#----------------------------------------------------------------------------- #
def get_fastq_input(wc):
    samps = config['samples'][wc.name]['fastq']

    if type(samps) is str:
        samps = [samps]

    infiles = [os.path.join(PATH_FASTQ, i) for i in samps]
    return(infiles)

def get_library_type(wc):
    lib = config['samples'][wc.name]['library']
    return(lib)

rule bowtie2:
    input:
        infile = get_fastq_input,
        genome = rules.bowtie2_build.output.outfile
    output:
        bamfile = os.path.join(PATH_MAPPED, "{name}/{name}.bam")
    params:
        threads  = 8,
        mem      ='16G',
        bowtie2  = SOFTWARE['bowtie2'],
        samtools = SOFTWARE['samtools'],
        library  = get_library_type,
        params_bowtie2 = PARAMS['bowtie2']
    log:
        log = os.path.join(PATH_LOG, "{name}.bowtie2.log")
    message:"""
        Mapping with bowtie2:
            sample: {input.infile}
            genome: {input.genome}
            output: {output.bamfile}
    """
    run:
        genome = input.genome.replace('.1.bt2','')
        if params.library in ['single','SINGLE']:
            map_args =  '-U ' + input.infile[0]

        if params.library in ['paired','PAIRED','pair','PAIR']:
            map_args = '-1 ' + input.infile[0] + ' -2 ' + input.infile[1]

        command = " ".join(
        [params.bowtie2,
        '-p', str(params.threads),
        '-x', genome,
        map_args,
        join_params("bowtie2", APP_PARAMS, params.params_bowtie2),
        '2>',log.log,
        '|', params.samtools,'view -bhS >', output.bamfile
        ])
        shell(command)

#----------------------------------------------------------------------------- #
rule samtools_sort:
    input:
        os.path.join(PATH_MAPPED, "{name}/{name}.bam")
    output:
        os.path.join(PATH_MAPPED, "{name}/{name}.sorted.bam")
    params:
        threads = 4,
        mem = '16G',
        samtools = SOFTWARE['samtools']
    message:"""
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
    message:"""
        Indexing bam file:\n
            input: {input}
    """
    shell:
        "{params.samtools} index {input}"


# ----------------------------------------------------------------------------- #
rule fastqc:
    input:
        infile = get_fastq_input
    output:
        outfile = os.path.join(PATH_QC, "{name}", "{name}.fastqc.done")
    params:
        outpath = os.path.join(PATH_QC, "{name}"),
        threads = 1,
        mem = '8G',
        fastqc = SOFTWARE['fastqc']
    message:"""
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

rule bam2bigWig:
    input:
        file = rules.bam2bed.output,
        chrlen = rules.chrlen.output
    output:
        file = os.path.join(os.getcwd(), PATH_MAPPED, "{name}/{name}.bw")
    params:
        threads  = 1,
        mem      = '16G',
        extend   = config['params']['extend'],
        scale    = config['params']['scale_bw']
    message:"""
        Making bigWig:
            input : {input.file}
            output: {output.file}
            scale:  {params.scale}
    """
    script:
        os.path.join(SCRIPT_PATH, 'BigWigExtend.R')

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
def get_files_macs(wc):
    paths = {}

    chip = config['peak_calling'][wc.name]['ChIP']
    if isinstance(chip,str):
        chip = [chip]
    chips = [os.path.join('Mapped','Bowtie', i, i + '.sorted.bam') for i in chip]
    paths['ChIP'] = chips
    
    cont = config['peak_calling'][wc.name]['Cont']
    if not cont == None:
        if isinstance(cont,str):
            cont = [cont]
        cont = [os.path.join('Mapped','Bowtie', i, i + '.sorted.bam') for i in cont]
        paths['Cont'] = cont

    return(paths)


rule macs2:
    input:
        unpack(get_files_macs)
    output:
        outfile = os.path.join(PATH_PEAK, "{name}", "{name}_peaks.narrowPeak")
    params:
        outpath = os.path.join(PATH_PEAK, "{name}"),
        name = "{name}",
        threads = 1,
        mem = '16G',
        macs2 = SOFTWARE['macs2'],
        params_macs = PARAMS['macs2']
    log:
        log = os.path.join(PATH_LOG, '{name}.macs.log')
    message:"""
        Running macs2:
            sample: {params.name}
            output: {output.outfile}
    """
    run:
        params_macs = params.params_macs
        if 'params' in config['peak_calling'][params.name].keys():
            if 'macs2' in config['peak_calling'][params.name]['params'].keys():
                params_macs.update(config['peak_calling'][params.name]['params']['macs2'])

        # checks whether the control samples are specified
        samples = ''
        samples = samples + " ".join(['-t'] + input.ChIP)
        
        if hasattr(input, 'Cont'):
            samples = samples + " ".join([' -c'] + input.Cont)
            
        command = " ".join(
        [params.macs2, 'callpeak',
        samples,
        '--outdir', params.outpath,
        '-n', params.name,
        join_params("macs2", APP_PARAMS, params_macs),
        '2>', log.log
        ])
        shell(command)

# ----------------------------------------------------------------------------- #
rule sort_narrow_peak:
    input:
        rules.macs2.output
    output:
        outfile = os.path.join(PATH_PEAK, "{name}", "{name}_qsort.narrowPeak")
    params:
        threads = 1,
        mem = '8G'
    message:"""
            Sorting narrow peak:
                input : {input}
                output: {output}
        """
    shell:"""
        sort -r -k9 -n {input} > {output.outfile}
    """
    
# ----------------------------------------------------------------------------- #
def get_chrlen(wc):
    matches = []
    for root, dirnames, filenames in os.walk(PATH_MAPPED):
        for filename in fnmatch.filter(filenames, 'chrlen.txt'):
            matches.append(os.path.join(root, filename))
    return(matches[0])
    
rule bedTobigBed:
    input:
        peaks  = rules.macs2.output,
        chrlen = get_chrlen
    output:
        outfile = os.path.join(PATH_PEAK, "{name}", "{name}.bb")
    params:
        threads = 1,
        mem = '8G',
        bedToBigBed = SOFTWARE['bedToBigBed']
    message:"""
            bedToBigBed:
                input : {input}
                output: {output}
        """
    shell:"""
        {params.bedToBigBed} -type=bed3+7 {input.peaks} {input.chrlen} {output.outfile} 
    """
        

# ----------------------------------------------------------------------------- #

def get_sample_idr(wc):
    name = config['idr'][wc.name]
    samps = dict(zip(name.keys(),[os.path.join(PATH_PEAK, i, i + '_qsort.narrowPeak') for i in name.values()]))
    return(samps)


rule idr:
    input:
        unpack(get_sample_idr)
    output:
        outfile = os.path.join(PATH_IDR, "{name}", "{name}.narrowPeak")
    params:
        threads = 1,
        mem = '8G',
        idr = SOFTWARE['idr'],
        params_idr = PARAMS['idr']
    log:
        log = os.path.join(PATH_LOG, '{name}.idr.log')
    message:"""
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
        '-l', log.log,
        '--plot',
        join_params("idr", APP_PARAMS, params.params_idr)
        ])
        shell(command)

# ----------------------------------------------------------------------------- #
if 'hub' in set(config.keys()):
    rule make_ucsc_hub:
        input:
            peaks  = BB,
            tracks = BW,
        output:
            outfile = os.path.join(PATH_HUB, HUB_NAME, 'done.txt')
        params:
            threads     = 1,
            mem         = '8G',
            hub         = config['hub'],
            genome_name = GENOME,
            paths       = TRACK_PATHS,
            path_hub    = os.path.join(PATH_HUB, HUB_NAME)
        log:
            log = os.path.join(PATH_LOG, 'UCSC_HUB.log')
        message:"""
                Running UCSC_HUB:
                    output: {output.outfile}
            """
        script:
            os.path.join(SCRIPT_PATH, 'Make_UCSC_HUB.R')
            

