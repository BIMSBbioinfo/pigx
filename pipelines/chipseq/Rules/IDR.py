# ----------------------------------------------------------------------------- #
def get_sample_idr(wc):
    name = config['idr'][wc.name]
    samps = dict(zip(name.keys(),[os.path.join(PATH_PEAK, i, i + '_qsort.narrowPeak') for i in name.values()]))
    return(samps)


rule idr:
    input:
        unpack(get_sample_idr)
    output:
        outfile = os.path.join(PATH_IDR, "{name}", "{name}.narrowPeak")
    params:
        threads = 1,
        mem = '8G',
        idr = SOFTWARE['idr'],
        params_idr = PARAMS['idr']
    log:
        log = os.path.join(PATH_LOG, '{name}.idr.log')
    message:"""
            Running IDR2:
                input : {input.ChIP1} {input.ChIP2}
                output: {output.outfile}
        """
    run:
        command = " ".join(
        [params.idr,
        '--samples', input.ChIP1, input.ChIP2,
        '--input-file-type', 'narrowPeak',
        '--rank', 'q.value',
        '--output-file', output.outfile,
        '-l', log.log,
        '--plot',
        join_params("idr", APP_PARAMS, params.params_idr)
        ])
        shell(command)