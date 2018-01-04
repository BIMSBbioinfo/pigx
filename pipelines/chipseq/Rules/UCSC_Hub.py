#----------------------------------------------------------------------------- #
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
def bedToBigBed_input(wc):
    suffix = get_macs2_suffix(wc.name, config)
    infile = os.path.join(PATH_PEAK, wc.name, wc.name + "_peaks." + suffix)
    return(infile)

rule bedTobigBed:
    input:
        peaks  = bedToBigBed_input,
        chrlen = rules.index_to_chrlen.output.outfile,
    output:
        outfile = os.path.join(PATH_PEAK, "{name}", "{name}.bb")
    params:
        threads     = 1,
        mem         = '8G',
        name        = '{name}',
        bedToBigBed = SOFTWARE['bedToBigBed']
    message:"""
            bedToBigBed:
                input : {input}
                output: {output}
        """
    run:
        suffix = get_macs2_suffix(params.name, config)
        if suffix == 'narrowPeak':
            ncol = 7
        else:
            ncol = 6

        command = ' '.join([
        params.bedToBigBed,
        '-type=bed3+' + str(ncol),
        input.peaks,
        input.chrlen,
        output.outfile
        ])
        shell(command)
