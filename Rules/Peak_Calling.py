# ----------------------------------------------------------------------------- #
def get_files_macs(wc):
    paths = {}

    chip = config['peak_calling'][wc.name]['ChIP']
    if isinstance(chip,str):
        chip = [chip]
    chips = [os.path.join('Mapped','Bowtie', i, i + '.sorted.bam') for i in chip]
    paths['ChIP'] = chips

    cont = config['peak_calling'][wc.name]['Cont']
    if not cont == None:
        if isinstance(cont,str):
            cont = [cont]
        cont = [os.path.join('Mapped','Bowtie', i, i + '.sorted.bam') for i in cont]
        paths['Cont'] = cont

    return(paths)


rule macs2:
    input:
        unpack(get_files_macs)
    output:
        outfile = os.path.join(PATH_PEAK, "{name}", "{name}_peaks.{class, narrow|broad}Peak")
    params:
        outpath = os.path.join(PATH_PEAK, "{name}"),
        name = "{name}",
        threads = 1,
        mem = '16G',
        macs2 = SOFTWARE['macs2'],
        params_macs = PARAMS['macs2']
    log:
        log = os.path.join(PATH_LOG, '{name}.macs.log')
    message:"""
        Running macs2:
            sample: {params.name}
            output: {output.outfile}
    """
    run:
        params_macs = params.params_macs
        if 'params' in config['peak_calling'][params.name].keys():
            if 'macs2' in config['peak_calling'][params.name]['params'].keys():
                params_macs.update(config['peak_calling'][params.name]['params']['macs2'])

        # checks whether the control samples are specified
        samples = ''
        samples = samples + " ".join(['-t'] + input.ChIP)

        if hasattr(input, 'Cont'):
            samples = samples + " ".join([' -c'] + input.Cont)

        command = " ".join(
        [params.macs2, 'callpeak',
        samples,
        '--outdir', params.outpath,
        '-n', params.name,
        join_params("macs2", APP_PARAMS, params_macs),
        '2>', log.log
        ])
        shell(command)

# ----------------------------------------------------------------------------- #
rule sort_peak:
    input:
        rules.macs2.output
    output:
        outfile = os.path.join(PATH_PEAK, "{name}", "{name}_qsort.{class, narrow|broad}Peak")
    params:
        threads = 1,
        mem = '8G'
    message:"""
            Sorting narrow peak:
                input : {input}
                output: {output}
        """
    shell:"""
        sort -r -k9 -n {input} > {output.outfile}
    """
