# ---------------------------------------------------------------------------- #
rule bowtie2_build:
    input:
        genome = GENOME_FASTA
    output:
        outfile = PREFIX + '.1.bt2'
    params:
        prefix  = PREFIX,
        threads = 1,
        mem = '32G',
        bowtie2_build = SOFTWARE['bowtie2-build']['executable']
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

# ---------------------------------------------------------------------------- #
rule index_to_chrlen:
        input:
            rules.bowtie2_build.output
        output:
            outfile = PREFIX + '.chrlen.txt'
        params:
            prefix  = PREFIX,
            threads = 1,
            mem     = '1G',
            bowtie2_inspect = SOFTWARE['bowtie2-inspect']['executable']
        log:
            os.path.join(PATH_LOG, "index_to_chrlen.log")
        message:
            """
                Constructing bowtie2 index:
                    input : {params.prefix}
                    output: {output.outfile}
            """
        shell:"""
            {params.bowtie2_inspect} -s {params.prefix} | grep Sequence | cut -f2,3 > {output.outfile} 2> {log}
        """

#----------------------------------------------------------------------------- #
rule bowtie2:
    input:
        infile = get_fastq_input,
        genome = rules.bowtie2_build.output.outfile
    output:
        bamfile = os.path.join(PATH_MAPPED, "{name}", "{name}.bam")
    params:
        threads        = 2,
        mem            =  '16G',
        bowtie2        = SOFTWARE['bowtie2']['executable'],
        samtools       = SOFTWARE['samtools']['executable'],
        library        = get_library_type,
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
        join_params("bowtie2", PARAMS, params.params_bowtie2),
        '2>',log.log,
        '|', params.samtools,'view -bhS >', output.bamfile
        ])
        shell(command)

#----------------------------------------------------------------------------- #
rule samtools_sort:
    input:
        os.path.join(PATH_MAPPED, "{name}", "{name}.bam")
    output:
        os.path.join(PATH_MAPPED, "{name}", "{name}.sorted.bam")
    params:
        threads = 4,
        mem = '16G',
        samtools = SOFTWARE['samtools']['executable']
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
        os.path.join(PATH_MAPPED, "{name}", "{name}.sorted.bam")
    output:
        os.path.join(PATH_MAPPED, "{name}", "{name}.sorted.bam.bai")
    params:
        threads = 1,
        mem = '8G',
        samtools = SOFTWARE['samtools']['executable']
    message:"""
        Indexing bam file:\n
            input: {input}
    """
    shell:
        "{params.samtools} index {input}"
