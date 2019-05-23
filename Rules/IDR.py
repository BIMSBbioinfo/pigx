# ----------------------------------------------------------------------------- #
def get_sample_idr(wc):
    name = config['idr'][wc.name]
    samps = [os.path.join(PATH_PEAK, i, i + '_qsort.bed') for i in name.values()]
    return(samps)


rule idr:
    input:
        unpack(get_sample_idr)
    output:
        outfile = os.path.join(PATH_IDR, "{name}", "{name}.bed")
    params:
        idr        = SOFTWARE['idr']['executable'],
        params_idr = PARAMS['idr']
    log:
        log = os.path.join(PATH_LOG, '{name}', '{name}.idr.log')
    message:"""
            Running IDR2:
                input : {input}
                output: {output.outfile}
        """
    run:
        command = " ".join(
        [params.idr,
        '--samples', input,
        '--input-file-type',  'bed',
        '--output-file-type', 'bed',
        '--rank', '9',
        '--output-file', output.outfile,
        '-l', log.log,
        '--plot',
        join_params("idr", PARAMS, params.params_idr)
        ])
        shell(command)
