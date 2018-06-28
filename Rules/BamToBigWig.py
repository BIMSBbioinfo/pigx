# ----------------------------------------------------------------------------- #
rule bam2bigWig:
    input:
        file    = os.path.join(PATH_MAPPED, "{name}","{name}" + BAM_SUFFIX),
        chrlen  = rules.index_to_chrlen.output.outfile
    output:
        outfile = os.path.join(PATH_MAPPED, "{name}", "{name}.bw")
    params:
        extend   = PARAMS['export_bigwig']['extend'],
        scale    = PARAMS['export_bigwig']['scale_bw'],
        Rscript  = SOFTWARE['Rscript']['executable'],
        library  = lambda wc: get_library_type(wc.name)
    log:
        logfile = os.path.join(PATH_LOG, 'bam2bigWig.log')
    message:"""
        Making bigWig:
            input : {input.file}
            output: {output.outfile}
            scale:  {params.scale}
    """
    run:
        RunRscript(input, output, params, log.logfile, 'BigWigExtend.R')

# ----------------------------------------------------------------------------- #
rule makelinks:
    input:
        file = rules.bam2bigWig.output.outfile
    output:
        outfile = os.path.join(PATH_BW, "{name}.bw")
    run:
        os.symlink(input.file, output.outfile)
