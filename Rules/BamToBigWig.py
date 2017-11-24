# ----------------------------------------------------------------------------- #
rule bam2bed:
    input:
        file   = os.path.join(PATH_MAPPED, "{name}/{name}.sorted.bam"),
        chrlen = rules.index_to_chrlen.output.outfile
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
        file   = rules.bam2bed.output,
        chrlen = rules.index_to_chrlen.output.outfile
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
