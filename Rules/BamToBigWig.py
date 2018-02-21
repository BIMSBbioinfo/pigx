# ----------------------------------------------------------------------------- #
rule bam2bed:
    input:
        file   = os.path.join(PATH_MAPPED, "{name}/{name}.sorted.bam"),
        chrlen = rules.index_to_chrlen.output.outfile
    output:
        outfile = os.path.join(PATH_MAPPED, "{name}", "{name}.bed")
    params:
        extend   = PARAMS['export_bigwig']['extend'],
        bamToBed = SOFTWARE['bamToBed']['executable']
    shell: """
        {params.bamToBed} -i {input.file} > {output.outfile}
    """

rule bam2bigWig:
    input:
        file    = rules.bam2bed.output.outfile,
        chrlen  = rules.index_to_chrlen.output.outfile
    output:
        outfile = os.path.join(PATH_MAPPED, "{name}", "{name}.bw")
    params:
        extend   = PARAMS['export_bigwig']['extend'],
        scale    = PARAMS['export_bigwig']['scale_bw'],
        Rscript  = SOFTWARE['Rscript']['executable']
    message:"""
        Making bigWig:
            input : {input.file}
            output: {output.outfile}
            scale:  {params.scale}
    """
    run:
        RunRscript(input, output, params, 'BigWigExtend.R')

# ----------------------------------------------------------------------------- #
rule makelinks:
    input:
        file = os.path.join(PATH_MAPPED, "{name}", "{name}" + '.bw')
    output:
        os.path.join(PATH_BW, "{name}.bw")
    shell: """
        ln -s {input.file} {output}
    """
