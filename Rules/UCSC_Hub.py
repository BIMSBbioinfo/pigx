#----------------------------------------------------------------------------- #
rule make_ucsc_hub:
        input:
            peaks  = BB,
            tracks = BW,
        output:
            outfile = os.path.join(PATH_HUB, HUB_NAME, 'done.txt')
        params:
            hub         = config['hub'],
            genome_name = GENOME,
            paths       = TRACK_PATHS,
            path_hub    = os.path.join(PATH_HUB, HUB_NAME),
            Rscript     = SOFTWARE['Rscript']['executable']
        log:
            logfile = os.path.join(PATH_LOG, 'UCSC_HUB.log')
        message:"""
                Running: UCSC_HUB:
                    output: {output.outfile}
            """
        run:
            RunRscript(input, output, params, log.logfile, 'Make_UCSC_HUB.R')


# ----------------------------------------------------------------------------- #
def bedToBigBed_input(wc):
    suffix = get_macs2_suffix(wc.name, CUSTOM_PARAMS)
    infile = os.path.join(PATH_PEAK, wc.name, wc.name + "_peaks." + suffix)
    return(infile)

rule bedTobigBed:
    input:
        peaks  = bedToBigBed_input,
        chrlen = GENOME_HASH[GENOME_TYPES['Main']]['genome_prefix'] + '.chrlen.txt'
    output:
        outfile = os.path.join(PATH_PEAK, "{name}", "{name}.bb")
    params:
        name        = '{name}',
        bedToBigBed = SOFTWARE['bedToBigBed']['executable']
    log:
        log = os.path.join(PATH_LOG, "{name}.bedtobigbed.log")
    message:"""
            bedToBigBed:
                input : {input}
                output: {output}
        """
    run:
        suffix = get_macs2_suffix(params.name, CUSTOM_PARAMS)
        if suffix == 'narrowPeak':
            ncol = 7
        else:
            ncol = 6

        command = ' '.join([
        params.bedToBigBed,
        '-type=bed3+' + str(ncol),
        input.peaks,
        input.chrlen,
        output.outfile,
        '2>',log.log])
        shell(command)
