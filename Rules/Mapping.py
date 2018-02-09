# ---------------------------------------------------------------------------- #
rule link_genome:
    input:
        genome  = GENOME_ORIG
    output:
        outfile = GENOME_FASTA
    params:
        prefix  = PREFIX,
        threads = 1,
        mem     = '1G',
    log:
        os.path.join(PATH_LOG, "link_genome.log")
    message:
        """
            Linking genome fasta:
                input : {input.genome}
                output: {output.outfile}
        """
    run:
        import os
        os.symlink(str(input.genome), str(output.outfile))
        
# ---------------------------------------------------------------------------- #
rule bowtie2_build:
    input:
        genome = GENOME_FASTA
    output:
        outfile = PREFIX + '.1.bt2'
    params:
        prefix  = PREFIX,
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
rule construct_genomic_windows:
        input:
            infile = rules.index_to_chrlen.output.outfile
        output:
            outfile = PREFIX + '.GenomicWindows.GRanges.rds'
        params:
            threads   = 1,
            mem       = '8',
            tilewidth = 10000,
            scriptdir = SCRIPT_PATH,
            Rscript   = SOFTWARE['Rscript']['executable']
        log:
            log = os.path.join(PATH_LOG, 'construct_genomic_windows.log')
        message:"""
                Running: construct_genomic_windows:
                    output: {output.outfile}
            """
        run:
            RunRscript(input, output, params, 'ConstructGenomicWindows.R')
            
#----------------------------------------------------------------------------- #
rule extract_nucleotide_frequency:
        input:
            genome_fasta    = GENOME_FASTA, 
            tilling_windows = rules.construct_genomic_windows.output.outfile
        output:
            outfile = PREFIX + '.NucleotideFrequency.GRanges.rds'
        params:
            threads   = 1,
            mem       = '8',
            scriptdir = SCRIPT_PATH,
            Rscript   = SOFTWARE['Rscript']['executable']
        log:
            log = os.path.join(PATH_LOG, 'extract_nucleotide_frequency.log')
        message:"""
                Running: extract_nucleotide_frequency:
                    output: {output.outfile}
            """
        run:
            RunRscript(input, output, params, 'Extract_Nucleotide_Frequency.R')

#----------------------------------------------------------------------------- #
rule bowtie2:
    input:
        infile = get_fastq_input,
        genome = rules.bowtie2_build.output.outfile
    output:
        bamfile = os.path.join(PATH_MAPPED, "{name}", "{name}.bam")
    params:
        threads        = config['execution']['rules']['bowtie2']['threads'],
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
        '-p', '{params.threads}',
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
        threads  = config['execution']['rules']['samtools_sort']['threads'],
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
        samtools = SOFTWARE['samtools']['executable']
    message:"""
        Indexing bam file:\n
            input: {input}
    """
    shell:
        "{params.samtools} index {input}"
