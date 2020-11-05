# ---------------------------------------------------------------------------- #
# Creates a symbolic link of the genome file into the project folder
# should be used inside rules
def link_genome(input, output):
    import os
    import magic as mg
    import gzip

    # ------------------------------------------------#
    # prevents errors during linking
    def trylink(infile, outfile):
        if os.path.exists(outfile):
            os.remove(outfile)
        try:
            os.symlink(infile, outfile)
        except:
            "Symbolic link exists"

    # ------------------------------------------------#
    # checks whether the input file is zipped, if zipped
    # unzips the file
    infile = str(input.genome)
    outfile = str(output.outfile)

    # this is necessary because the [g]unzip command
    # dies if the existing file is in the directory
    # necessary if the reference gets updated
    if os.path.exists(outfile):
        os.remove(outfile)

    genome_file_type = mg.from_file(infile, mime=True)

    # this is ugly piece of code
    if genome_file_type.find('gzip') > 0:

        with gzip.open(infile, 'r') as file:
            with open(outfile, 'w') as out:
                for line in file:
                    out.write(str(line.decode('utf-8')))

    elif genome_file_type == 'text/plain':
        trylink(infile, outfile)

    else:
        os.sys.exit('Unknow input genome file format')

# ---------------------------------------------------------------------------- #
rule link_genome:
    input:
        genome  = lambda wc: GENOME_HASH[wc.genome_type]['genome_location']
    output:
        outfile = os.path.join(PATH_INDEX, '{genome_type}','{genome}.fa')
    wildcard_constraints:
        genome_type = GENOME_TYPES_CONSTRAINT
    params:
        threads = config['execution']['rules']['link_genome']['threads'],
        mem     = config['execution']['rules']['link_genome']['memory'],
    log:
        logfile = os.path.join(PATH_LOG, "link_{genome_type}_genome_{genome}.log")
    message:
        """
            Linking genome fasta:
                input : {input.genome}
                output: {output.outfile}
        """
    run:
        link_genome(input, output)


# ---------------------------------------------------------------------------- #
rule bowtie2_build:
    input:
        genome  = lambda wc: GENOME_HASH[wc.genome_type]['genome_link']
    output:
        outfile = os.path.join(PATH_INDEX, '{genome_type}','{genome}.1.bt2')
    wildcard_constraints:
        genome_type = GENOME_TYPES_CONSTRAINT
    params:
        bowtie2_build = SOFTWARE['bowtie2-build']['executable']
    log:
        log = os.path.join(PATH_LOG, "bowtie2_build_{genome_type}_genome_{genome}.log")
    message:
        """
            Constructing bowtie2 index:
                input : {input.genome}
                output: {output.outfile}
        """
    run:
        prefix = output.outfile.replace('.1.bt2','')

        command = " ".join(
        [params.bowtie2_build,
        input.genome,
        prefix,
        '>>', log.log,
        '2>&1'
        ])
        shell(command)

# ---------------------------------------------------------------------------- #
rule index_to_chrlen:
        input:
            infile = rules.bowtie2_build.output
        output:
            outfile = os.path.join(PATH_INDEX, '{genome_type}','{genome}.chrlen.txt')
        wildcard_constraints:
            genome_type = GENOME_TYPES_CONSTRAINT
        params:
            bowtie2_inspect = SOFTWARE['bowtie2-inspect']['executable'],
            grep = SOFTWARE['grep']['executable'],
            cut = SOFTWARE['cut']['executable']
        log:
            log = os.path.join(PATH_LOG, "index_to_chrlen_{genome_type}_genome_{genome}.log")
        message:
            """
                Extracting chromosome lengths from index:
                    input : {input.infile}
                    output: {output.outfile}
            """
        run:
            prefix = str(input.infile).replace('.1.bt2','')

            command = " ".join(
            [params.bowtie2_inspect,
            '-s', prefix,
            '|',  params.grep, 'Sequence',
            '|',  params.cut, '-f2,3',
            '>',  '{output.outfile}',
            '2>', log.log
            ])
            shell(command)


#----------------------------------------------------------------------------- #
rule construct_genomic_windows:
        input:
            infile = rules.index_to_chrlen.output.outfile
        output:
            outfile = os.path.join(PATH_INDEX, "{genome_type}",  '{genome}.GenomicWindows.GRanges.rds')
        wildcard_constraints:
            genome_type = GENOME_TYPES_CONSTRAINT
        params:
            threads   = 1,
            mem       = '8',
            tilewidth = 10000,
            scriptdir = SCRIPT_PATH,
            Rscript   = SOFTWARE['Rscript']['executable']
        log:
            logfile = os.path.join(PATH_LOG, 'construct_genomic_windows_{genome_type}_genome_{genome}.log')
        message:"""
                Running: construct_genomic_windows:
                    output: {output.outfile}
            """
        run:
            RunRscript(input, output, params, log.logfile, 'ConstructGenomicWindows.R')

#----------------------------------------------------------------------------- #
rule extract_nucleotide_frequency:
        input:
            genome_fasta    = GENOME_FASTA,
            tilling_windows = rules.construct_genomic_windows.output.outfile
        output:
            outfile = os.path.join(PATH_INDEX, "{genome_type}", '{genome}.NucleotideFrequency.GRanges.rds')
        wildcard_constraints:
            genome_type = GENOME_TYPES_CONSTRAINT
        params:
            threads   = 1,
            mem       = '8',
            scriptdir = SCRIPT_PATH,
            Rscript   = SOFTWARE['Rscript']['executable']
        log:
            logfile = os.path.join(PATH_LOG, 'extract_nucleotide_frequency_{genome_type}_genome_{genome}.log')
        message:"""
                Running: extract_nucleotide_frequency:
                    output: {output.outfile}
            """
        run:
            RunRscript(input, output, params, log.logfile, 'Extract_Nucleotide_Frequency.R')

#----------------------------------------------------------------------------- #
def get_trimmed_files(sample):
    infiles = TRIM_GALORE_FILES[sample]
    return(infiles)

rule bowtie2:
    input:
        infile = lambda wc: get_trimmed_files(wc.name),
        genome = lambda wc: GENOME_HASH[wc.genome_type]['bowtie_index']
    output:
        bamfile = os.path.join(PATH_MAPPED, "{genome_type}", "{name}", "{name}.bam")
    wildcard_constraints:
        genome_type = GENOME_TYPES_CONSTRAINT
    params:
        threads        = config['execution']['rules']['bowtie2']['threads'],
        bowtie2        = SOFTWARE['bowtie2']['executable'],
        samtools       = SOFTWARE['samtools']['executable'],
        library        = lambda wc: get_library_type(wc.name),
        params_bowtie2 = PARAMS['bowtie2']
    log:
        logfile = os.path.join(PATH_LOG, "{name}", "{name}.bowtie2_{genome_type}.log")
    message:"""
        Mapping with bowtie2:
            sample: {input.infile}
            genome: {input.genome}
            output: {output.bamfile}
    """
    run:
        genome = input.genome.replace('.1.bt2','')
        if params.library in ['single','SINGLE']:
            map_args = '-U {}'.format(','.join(input.infile))

        if params.library in ['paired','PAIRED','pair','PAIR']:
            map_args = '-1 {} -2 {}'.format(','.join(input.infile[::2]),','.join(input.infile[1::2]))

        command = " ".join(
        [params.bowtie2,
        '-p', '{params.threads}',
        '-x', genome,
        map_args,
        join_params("bowtie2", PARAMS, params.params_bowtie2),
        '2>',log.logfile,
        '|', params.samtools,'view -bhS >', output.bamfile
        ])
        shell(command)

#----------------------------------------------------------------------------- #
# TODO: once samtools 1.10 is available via guix:
#   -   optionally tag optical duplicates '-d 100'
#   -   save stats to file with '-f'
rule samtools_markduplicate:
    input:
        os.path.join(PATH_MAPPED, "{genome_type}", "{name}", "{name}.bam")
    output:
        os.path.join(PATH_MAPPED, "{genome_type}", "{name}", "{name}.sorted.bam")
    params:
        threads  = config['execution']['rules']['samtools_sort']['threads'],
        samtools = SOFTWARE['samtools']['executable']
    log:
        logfile = os.path.join(PATH_LOG, "{name}","{name}.{genome_type}_genome.samtools_markdup.log")
    message:"""
            Marking duplicate reads:
                input: {input}
                output: {output}
        """
    shell: """
        {params.samtools} sort --threads {params.threads} -n {input} 2> {log}| \
        {params.samtools} fixmate --threads {params.threads}  -m - - 2> {log}| \
        {params.samtools} sort --threads {params.threads}  - 2> {log}| \
        {params.samtools} markdup --threads {params.threads} -S -s - {output} \
        2> {log}
    """


# ----------------------------------------------------------------------------- #
rule samtools_index:
    input:
        os.path.join(PATH_MAPPED, "{genome_type}","{name}", "{name}.sorted.bam" )
    output:
        os.path.join(PATH_MAPPED, "{genome_type}","{name}", "{name}.sorted.bam.bai")
    wildcard_constraints:
        genome_type = GENOME_TYPES_CONSTRAINT
    params:
        samtools = SOFTWARE['samtools']['executable']
    log:
        logfile = os.path.join(PATH_LOG, "{name}","{name}.{genome_type}_genome.samtools_index.log")
    message:"""
        Indexing bam file:\n
            input: {input}
    """
    shell:
        "{params.samtools} index {input} 2> {log}"


rule samtools_idxstats:
  input:
    os.path.join(PATH_MAPPED, "{genome_type}","{name}", "{name}.sorted.bam")
  output:
    os.path.join(PATH_MAPPED, "{genome_type}","{name}", "{name}.sorted.idxstats.txt")
  wildcard_constraints:
    genome_type = GENOME_TYPES_CONSTRAINT
  params:
    samtools = SOFTWARE['samtools']['executable']
  log:
    logfile = os.path.join(PATH_LOG, "{name}","{name}.{genome_type}_genome.samtools_idxstats.log")
  shell:
        "{params.samtools} idxstats {input} > {output} 2> {log}"


rule samtools_stat:
  input:
    os.path.join(PATH_MAPPED, "{genome_type}","{name}", "{name}.sorted.bam")
  output:
    os.path.join(PATH_MAPPED, "{genome_type}","{name}", "{name}.sorted.stats.txt")
  wildcard_constraints:
    genome_type = GENOME_TYPES_CONSTRAINT
  params:
    samtools = SOFTWARE['samtools']['executable']
  log:
    logfile = os.path.join(PATH_LOG, "{name}","{name}.{genome_type}_genome.samtools_stats.log")
  shell:
        "{params.samtools} stats {input} > {output} 2> {log}"


rule samtools_flagstat:
  input:
    os.path.join(PATH_MAPPED, "{genome_type}","{name}", "{name}.sorted.bam")
  output:
    os.path.join(PATH_MAPPED, "{genome_type}","{name}", "{name}.sorted.flagstat.txt")
  wildcard_constraints:
    genome_type = GENOME_TYPES_CONSTRAINT
  params:
    samtools = SOFTWARE['samtools']['executable']
  log:
    logfile = os.path.join(PATH_LOG, "{name}","{name}.{genome_type}_genome.samtools_flagstat.log")
  shell:
        "{params.samtools} flagstat {input} > {output} 2> {log}"
