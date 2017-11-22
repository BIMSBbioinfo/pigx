"""
jump conda3
source ./activate p35
SNAKEFILE='/home/vfranke/Projects/AAkalin_PIX/Snake_ChIPseq.py'
WORKDIR='/data/akalin/vfranke/AAkalin_PIX/ChIP'
CONFIGFILE='/home/vfranke/Projects/AAkalin_PIX/config.yaml'
PATH='/home/vfranke/bin/Software/miniconda3/envs/p35/bin:/usr/local/bin:/usr/bin:/bin:/home/vfranke/.guix-profile/bin:/home/vfranke/.guix-profile/sbin:/home/vfranke/bin'
# Beast run

snakemake -R --snakefile $SNAKEFILE --directory $WORKDIR --jobs 4 --rerun-incomplete --configfile $CONFIGFILE --latency-wait 30 --dryrun

snakemake -R --snakefile $SNAKEFILE --directory $WORKDIR --jobs 4 --rerun-incomplete --configfile $CONFIGFILE --latency-wait 30

snakemake -R --snakefile $SNAKEFILE --directory $WORKDIR --jobs 4 --rerun-incomplete --configfile $CONFIGFILE --latency-wait 30 --dryrun

# Visualize the graph horizontally
snakemake --dag --snakefile $SNAKEFILE --directory $WORKDIR --configfile $CONFIGFILE | perl -pe 's|graph\[|graph[rankdir=LR,|'| dot -T png > Snake_ChIPseq.png

# Cluster run
snakemake -R --snakefile $SNAKEFILE --directory $WORKDIR --timestamp --jobs 24 --cluster "qsub -V -l h_vmem={params.mem} -pe smp {params.threads} -l h_rt=36:00:00" --rerun-incomplete --latency-wait 30

for i in `ls $WORKDIR | grep -v Fastq$ | grep -v Genome`;do rm -r $WORKDIR/$i;done

rm $(snakemake --snakefile $SNAKEFILE --directory $WORKDIR --configfile $CONFIGFILE --summary | tail -n+2 | cut -f1)


gp -i gcc-toolchain@5 gfortran@5 pkg-config curl cairo libxt libxml2 openssl xz zlib -p pix

gp -i gcc -p pix

bash
cd /home/vfranke/bin/guix
profile=/home/vfranke/bin/guix/pix
source $profile/etc/profile

export PATH="$profile/bin${PATH:+:}$PATH"
export PKG_CONFIG_PATH="pix/lib/pkgconfig${PKG_CONFIG_PATH:+:}$PKG_CONFIG_PATH"
export XDG_DATA_DIRS="pix/share${XDG_DATA_DIRS:+:}$XDG_DATA_DIRS"
export GIO_EXTRA_MODULES="pix/lib/gio/modules${GIO_EXTRA_MODULES:+:}$GIO_EXTRA_MODULES"
export R_LIBS_SITE="pix/site-library/${R_LIBS_SITE:+:}$R_LIBS_SITE"
export GUIX_LOCPATH=$profile/lib/locale
export CURL_CA_BUNDLE="/etc/ssl/certs/ca-bundle.crt"
/home/vfranke/bin/guix/pix/bin/R
source("https://bioconductor.org/biocLite.R")
biocLite('S4Vectors')
biocLite(c('data.table','genomation','GenomicRanges','RCAS','rtracklayr','stringr', 'plotly', 'DT', 'data.table', 'topGO', 'motifRG', 'biomaRt', 'AnnotationDbi', 'GenomicRanges', 'BSgenome.Hsapiens.UCSC.hg19', 'GenomeInfoDb', 'Biostrings', 'rtracklayer', 'org.Hs.eg.db', 'GenomicFeatures', 'genomation', 'rmarkdown', 'knitr', 'S4Vectors'), dependencies=TRUE)




Integrate_Expression_Mouse_Hamster_Difference
1. Figure out how to prevent bombing

### 13. config file check

Q:
1. How to write dynamic rules? - that the rules execute based on the input parameters?
    - if the user supplies the genome index, do not make the index
2. What to do downstream? Some basic plots? ChromHMM
3. How to speed up the extension?

"""

SCRIPT_PATH = os.path.join(workflow.basedir,'Scripts')
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

# Software executables
SOFTWARE     = SOFTWARE_CONFIG['software']
# ---------------------------------------------------------------------------- #
# Variable definition
GENOME       = config['genome']['name']
GENOME_FASTA = config['genome']['fasta']
NAMES        = config['samples'].keys()
PATH_FASTQ   = config['fastq']
PARAMS       = config['params']
ANNOTATION   = config['annotation']


# Directory structure definition
PATH_ANNOTATION = "Annotation"
PATH_MAPPED     = "Mapped/Bowtie"
PATH_QC         = "FastQC"
PATH_INDEX      = 'Bowtie2_Index'
PATH_LOG        = 'Log'
PATH_PEAK       = 'Peaks/MACS2'
PATH_BW         = 'BigWig'
PATH_IDR        = 'IDR'
PATH_HUB        = 'UCSC_HUB'
PATH_FEATURE    = "Analysis/Feature_Combination"

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
COMMAND    = []
INDEX      = [PREFIX + '.1.bt2']
BOWTIE2    = expand(os.path.join(PATH_MAPPED, "{name}", "{name}.sorted.bam.bai"), name=NAMES)
FASTQC     = expand(os.path.join(PATH_QC,     "{name}", "{name}.fastqc.done"), name=NAMES)
BW         = expand(os.path.join(os.getcwd(), PATH_MAPPED, "{name}", "{name}.bw"),  name=NAMES)
LINKS      = expand(os.path.join(PATH_BW,  "{ex_name}.bw"),  ex_name=NAMES)


COMMAND = COMMAND + INDEX + BOWTIE2 + BW + LINKS + FASTQC 
# ----------------------------------------------------------------------------- #
if 'peak_calling' in set(config.keys()):
    MACS  = []
    QSORT = []
    suffix = 'narrowPeak'
    for name in config['peak_calling'].keys():
        suffix = get_macs2_suffix(name, config)

        MACS    = MACS  + [os.path.join(PATH_PEAK,  name, name + "_peaks." + suffix)]
        QSORT   = QSORT + [os.path.join(PATH_PEAK,  name, name + "_qsort." + suffix)]

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

# ----------------------------------------------------------------------------- #
if 'feature_combination' in set(config.keys()):
    peak_files = []
    if 'idr' in config['feature_combination'].keys():
        peak_files = peak_files + [os.path.join(PATH_IDR,  x, x + '.narrowPeak') for x in config['feature_combination']['idr']]

    if 'peaks' in config['feature_combination'].keys():
        peak_files = peak_files + [os.path.join(PATH_PEAK, x, x + '_peaks.narrowPeak') for x in config['feature_combination']['peaks']]
    FEATURE = [os.path.join(PATH_FEATURE,'Feature_Combination.tsv')]
    COMMAND = COMMAND + FEATURE

if 'gtf' in set(config['annotation'].keys)
    PREPARE_ANNOTATION = [os.path.join(PATH_ANNOTATION,'Processed_Annotation.rds')]
    EXTRACT_SIGNAL     = expand(os.path.join(PATH_PEAK,'{name}','{name}.Extract_Signal.rds'), name=NAMES)
    ANNOTATE_PEAKS     = expand(os.path.join(PATH_PEAK,'{name}','{name}.Annotate_Peaks.rds'), name=NAMES)
        
    COMMAND = COMMAND + PREPARE_ANNOTATION + EXTRACT_SIGNAL + ANNOTATE_PEAKS

# ----------------------------------------------------------------------------- #
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
        threads  = 2,
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
        outfile = os.path.join(PATH_PEAK, "{name}", "{name}_peaks.{class, narrow|broad}Peak")
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
rule sort_peak:
    input:
        rules.macs2.output
    output:
        outfile = os.path.join(PATH_PEAK, "{name}", "{name}_qsort.{class, narrow|broad}Peak")
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

def bedToBigBed_input(wc):
    suffix = get_macs2_suffix(wc.name, config)
    infile = os.path.join(PATH_PEAK, wc.name, wc.name + "_peaks." + suffix)
    return(infile)

rule bedTobigBed:
    input:
        peaks  = bedToBigBed_input,
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
                Running: UCSC_HUB:
                    output: {output.outfile}
            """
        script:
            os.path.join(SCRIPT_PATH, 'Make_UCSC_HUB.R')

# ----------------------------------------------------------------------------- #
if 'feature_combination' in set(config.keys()):
    rule feature_combination:
        input:
            features = peak_files,
            bw = BW
        output:
            outfile = os.path.join(PATH_FEATURE,'Feature_Combination.tsv')
        params:
            threads     = 1,
            mem         = '8G',
            annotation  = ANNOTATION,
            outpath     = PATH_FEATURE,
            scriptdir   = SCRIPT_PATH
        log:
            log = os.path.join(PATH_LOG, 'feature_combination.log')
        message:"""
                Running: feature_combination:
                    output: {output.outfile}
            """
        script:
            os.path.join(SCRIPT_PATH, 'Feature_Combinaton.R')


#----------------------------------------------------------------------------- #
# Signal processing for the downstream analysis
rule prepare_annotation:
        input:
            annotation = ANNOTATION,
        output:
            outfile = os.path.join(PATH_ANNOTATION,'Processed_Annotation.rds')
        params:
            threads     = 1,
            mem         = '16G',
            scriptdir = SCRIPT_PATH
        log:
            log = os.path.join(PATH_LOG, 'prepare_annotation.log')
        message:"""
                Running: prepare_annotation:
                    output: {output.outfile}
            """
        script:
            os.path.join(SCRIPT_PATH, 'Prepare_Annotation.R')
            
            
#----------------------------------------------------------------------------- #
rule extract_signal:
        input:
            annotation = ANNOTATION,
            bed        = rules.sort_peak.outfile,
            wig        = rules.bam2bigWig
        output:
            outfile    = os.path.join(PATH_PEAK,'{name}','{name}.Extract_Signal.rds')
        params:
            threads     = 1,
            mem         = '16G',
            expand_peak = config['params']['rule extract_signal']['expand_peak']
            bin_num     = config['params']['rule extract_signal']['bin_num']
            scriptdir   = SCRIPT_PATH
        log:
            log = os.path.join(PATH_LOG, 'prepare_annotation.log')
        message:"""
                Running: extract_signal:
                    bed: {input:bed}
                    wig: {input:wig}
                    output: {output.outfile}
            """
        script:
            os.path.join(SCRIPT_PATH, 'Extract_Signal.R') 
            
#----------------------------------------------------------------------------- #
# Annotate Peaks
rule annotate_peaks:
        input:
            annotation = ANNOTATION,
            peaks      = rules.sort_peak.outfile,
        output:
            outfile    = os.path.join(PATH_PEAK,'{name}','{name}.Annotate_Peaks.rds')
        params:
            threads     = 1,
            mem         = '16G',
            expand_peak = config['params']['rule extract_signal']['expand_peak'],
            bin_num     = config['params']['rule extract_signal']['bin_num'],
            peakname    = '{name}',
            scriptdir   = SCRIPT_PATH
        log:
            log = os.path.join(PATH_LOG, 'prepare_annotation.log')
        message:"""
                Running: annotate_peaks:
                    bed: {input:bed}
                    wig: {input:wig}
                    output: {output.outfile}
            """
        script:
            os.path.join(SCRIPT_PATH, 'Annotate_Peaks.R')
