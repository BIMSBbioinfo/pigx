"""
SNAKEFILE='/home/vfranke/Projects/AAkalin_PIX/Scripts/ChIP/Snake_ChIPseq.py'
WORKDIR='/data/akalin/vfranke/AAkalin_PIX/ChIP'
CONFIGFILE='/home/vfranke/Projects/AAkalin_PIX/Scripts/ChIP/config.yaml'

# Beast run
snakemake -R --snakefile $SNAKEFILE --directory $WORKDIR --jobs 4 --rerun-incomplete --configfile $CONFIGFILE --latency-wait 30

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
12. Add sampleshee to yaml file
### 13. config file check

Q:
1. How to write dynamic rules? - that the rules execute based on the input parameters?
    - if the user supplies the genome index, do not make the index
2. What to do downstream? Some basic plots? ChromHMM
3. How to speed up the extension?

"""
import glob
import os
import re
from snakemake.utils import R

localrules: makelinks

# Variable definition
GENOME = config['genome']
NAMES = config['samples'].keys()

# Directory structure definition
PATH_MAPPED = "Mapped/Bowtie"
PATH_QC     = "FastQC"
PATH_GENOME = 'Bowtie_Index'
PATH_LOG    = 'Log'
PATH_PEAK   = 'Peaks/MACS2'
PATH_BW     = 'Links_bw'


# ---------------------------------------------------------------------------- #
# config defaults
if not ('extend' in config['params'].keys()):
    config['params']['extend'] = 0

PREFIX=''
if 'index' in config.keys() and config['index'] != None:
    PREFIX = config['index']
else:
    PREFIX = os.path.join(PATH_GENOME, os.path.splitext(os.path.basename(GENOME))[0])

# ---------------------------------------------------------------------------- #
# config check
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def check_config(config):

    message = ''

    # checks for proper top level config categories
    params = ['genome','index','params','samples','peak_calling']
    params_diff = set(config.keys()) - set(params)
    if len(params_diff) > 0:
        message = message + "config file contains dissalowed categories\n"

    # checks for correspondence between peak calling and samples
    samples = list(config['samples'].keys())
    peaks = [config['peak_calling'][i].values() for i in list(config['peak_calling'].keys())]
    samples_diff = (set(peaks[0]) - set(samples))
    if len(samples_diff) > 0:
        message = message + "some peak calling samples are not specified\n"

    # checks for index or genome specification
    if len(config['genome']) == 0 and len(config['index']):
        message = message + "neither genome nor index are specified\n"

    # checks whether extend is a number
    # if not (is.number(config['params']['extend'])):
    #     message = message + "extend must be a number\n"

    if len(message) > 0:
        print(message)
        return(1);

    return(0)

if check_config(config) == 1:
    print('config.yaml is not properly formatted - exiting')
    quit()



# ----------------------------------------------------------------------------- #
# RULE ALL
INDEX   = [PREFIX + '.1.ebwt']
MAPPING = expand(os.path.join(PATH_MAPPED, "{name}", "{name}.sorted.bam.bai"), name=NAMES)
FASTQC  = expand(os.path.join(PATH_QC,     "{name}", "{name}.fq_fastqc.zip"), name=NAMES)
BW      = expand(os.path.join(os.getcwd(), PATH_MAPPED, "{name}", "{name}.bw"),  name=NAMES)
LINKS   = expand(os.path.join(PATH_BW,  "{ex_name}.bw"),  ex_name=NAMES)
MACS    = expand(os.path.join(PATH_PEAK,    "{name}", "{name}.narrowPeak"), name=config['peak_calling'].keys())

rule all:
    input:
        INDEX + MAPPING + BW + LINKS + FASTQC + MACS


# ----------------------------------------------------------------------------- #
rule bowtie_build:
    input:
        genome = GENOME
    output:
        outfile = PREFIX + '.1.ebwt'
    params:
        prefix = PREFIX,
        threads = 1,
        mem = '32G'
    log:
        os.path.join(PATH_LOG, "bowtie_build.log")
    message: """
        Bowtie build ...
    """
    shell:     """
        bowtie-build {input.genome} {params.prefix} 2> {log}
    """

#----------------------------------------------------------------------------- #
def fastq_input(wc):
    samps = config['samples'][wc.name]

    if type(samps) is str:
        samps = [samps]

    infiles = [os.path.join('Fastq', i) for i in samps]
    return(infiles)

rule bowtie:
    input:
        infile = fastq_input,
        genome = rules.bowtie_build.output.outfile
    output:
        bamfile = os.path.join(PATH_MAPPED, "{name}/{name}.bam")
    params:
        threads=8,
        mem='16G'
    log:
        log = os.path.join(PATH_MAPPED, "{name}/{name}.bowtie.log")
    message: """
        Mapping ...
    """
    run:
        genome = input.genome.replace('.1.ebwt','')

        command = input.infile[0] + ' | bowtie -p ' + str(params.threads) +' -S -k 1 -m 1 --best --strata ' + genome  + ' -  2> ' + log.log + ' | samtools view -bhS > ' + output.bamfile

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
        mem = '16G'
    message: """
        Sorting ...
    """
    shell: """
        samtools sort --threads {params.threads} -o {output} {input}
    """

# ----------------------------------------------------------------------------- #
rule samtools_index:
    input:
        os.path.join(PATH_MAPPED, "{name}/{name}.sorted.bam")
    output:
        os.path.join(PATH_MAPPED, "{name}/{name}.sorted.bam.bai")
    params:
        threads = 1,
        mem = '8G'
    message: """
        Index ...
    """
    shell:
        "samtools index {input}"


# ----------------------------------------------------------------------------- #
rule fastqc:
    input:
        infile = fastq_input
    output:
        outfile = os.path.join(PATH_QC, "{name}", "{name}.fq_fastqc.zip")
    params:
        outpath = os.path.join(PATH_QC, "{name}"),
        threads = 1,
        mem = '8G'
    message:
        """FastQC"""
    shell: """
        fastqc --outdir {params.outpath} --extract -f fastq -t 1 {input.infile}
    """

# ----------------------------------------------------------------------------- #
rule chrlen:
    input:
        os.path.join(PATH_MAPPED, "{name}/{name}.sorted.bam")
    output:
        os.path.join(PATH_MAPPED, "{name}/chrlen.txt")
    params:
        threads = 1,
        mem = '4G'
    message: """
        chrlen ...
    """
    shell:"""
        samtools view -H {input} | grep @SQ | perl -pe 's/^@.+?://;s/LN://' > {output}
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
        mem = '16G'
    message: """
        bam2bed extend...
    """
    shell: """
        bamToBed -i {input.file} > {output}
    """

# perl -lane 'if(@F[5] eq '+'){{print join     (\"\t\",@F[0..1],@F[1]+{params.extend},@F[3..5]);}}else{{print join (\"\t\",@F[0],@F[2]-{params.extend},@F[2],@F[3..5]);}}' > {output}

rule bam2bed_extend:
    input:
        file = rules.bam2bed.output,
        chrlen = rules.chrlen.output
    output:
        os.path.join(PATH_MAPPED, "{name}/{name}.ext.bed")
    params:
        extend=config['params']['extend'],
        threads = 1,
        mem = '16G'
    message: """
        bam2bed extend...
    """
    shell:"""
    R --silent -e "library(data.table); f = fread('{input.file}'); setnames(f,c('chr','start','end','name','x','strand')); f[strand == '+', end := start[strand == '+'] + as.integer({params.extend})]; f[strand == '-', start := end[strand == '-'] - as.integer({params.extend})]; write.table(f, '{output}', row.names=F, col.names=F, quote=F, sep='\\t');"
    """
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
        bedgraph = os.path.join(PATH_MAPPED, "{name}/{name}.bedGraph")
    message: """
        bam2bedgraph ...
    """
    shell: """
        genomeCoverageBed -i {input.bedfile} -g {input.chrlen} -bg > {params.bedgraph}
        wigToBigWig {params.bedgraph} {input.chrlen} {output}
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
def get_sample_name(wc):
    name = config['peak_calling'][wc.name]
    samps = dict(zip(name.keys(),[os.path.join('Mapped','Bowtie',i, i + '.sorted.bam') for i in name.values()]))
    return(samps)

rule macs:
    input:
        unpack(get_sample_name)
    output:
        outdir = os.path.join(PATH_PEAK, "{name}", "{name}.narrowPeak")
    params:
        name = "{name}",
        threads = 1,
        mem = '16G'
    log:
    	os.path.join(PATH_LOG, 'macs.log')
    message: """
        macs ...
    """
    shell: """
        macs2 callpeak -t {input.ChIP} -c {input.Cont} --outdir {output.outdir} -n {params.name} 2> {log}
    """


