# ----------------------------------------------------------------------------- #
rule bam2bigWig:
    input:
        # bam files
        bamfile = os.path.join(PATH_MAPPED, "{genome_type}", "{name}","{name}" + BAM_SUFFIX),
        # chromosome lengths
        chrlen  = lambda wc: GENOME_HASH[wc.genome_type]['genome_prefix'] + '.chrlen.txt',
        # bowtie mapping stats
        stats   = os.path.join(PATH_RDS, "BowtieLog.rds")
    output:
        outfile = os.path.join(PATH_MAPPED, "{genome_type}", "{name}", "{name}.bw")
    params:
        extend      = PARAMS['export_bigwig']['extend'],
        scale       = PARAMS['export_bigwig']['scale_bw'],
        Rscript     = SOFTWARE['Rscript']['executable'],
        library     = lambda wc: get_library_type(wc.name),
        spikein     = lambda wc: get_spikein_information(wc.name),
        sample_name = lambda wc: wc.name
    log:
        logfile = os.path.join(PATH_LOG, "{name}", 'bam2bigWig.log')
    message:"""
        Making bigWig:
            input :   {input.bamfile}
            output:   {output.outfile}
            scale:    {params.scale}
            spikein:  {params.spikein}
    """
    run:
        RunRscript(input, output, params, log.logfile, 'BigWigExtend.R')

# ----------------------------------------------------------------------------- #
rule makelinks:
    input:
        file = os.path.join(PATH_MAPPED, "{genome_type}", "{name}", "{name}.bw")
    output:
        outfile = os.path.join(PATH_BW, "{genome_type}", "{name}.bw")
    message:"""
        Linking bigWig:
            input :   {input.file}
            output:   {output.outfile}
    """
    run:
        trylink(input.file, output.outfile)
