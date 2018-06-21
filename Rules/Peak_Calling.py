# ----------------------------------------------------------------------------- #
def get_files_macs(wc):
    paths = {}

    chip = config['peak_calling'][wc.name]['ChIP']
    if isinstance(chip,str):
        chip = [chip]
    chips = [os.path.join(PATH_MAPPED, i, i + '.sorted.bam') for i in chip]
    paths['ChIP'] = chips

    cont = config['peak_calling'][wc.name]['Cont']
    if not cont == None:
        if isinstance(cont,str):
            cont = [cont]
        cont = [os.path.join(PATH_MAPPED, i, i + '.sorted.bam') for i in cont]
        paths['Cont'] = cont

    return(paths)


rule macs2:
    input:
        unpack(get_files_macs),
        chrlen = rules.index_to_chrlen.output
    output:
        outfile = os.path.join(PATH_PEAK, "{name}", "{name}_peaks.{class, narrow|broad}Peak")
    params:
        outpath     = os.path.join(PATH_PEAK, "{name}"),
        name        = "{name}",
        macs2       = SOFTWARE['macs2']['executable'],
        params_macs = PARAMS['macs2']
    log:
        logfile = os.path.join(PATH_LOG, '{name}.macs.log')
    message:"""
        Running macs2:
            sample: {params.name}
            output: {output.outfile}
    """
    run:
        params_macs = params.params_macs
        if params.name in set(CUSTOM_PARAMS.keys()):
            if not CUSTOM_PARAMS[params.name] == None:
                if 'macs2' in CUSTOM_PARAMS[params.name].keys():
                    params_macs.update(CUSTOM_PARAMS[params.name]['macs2'])

        # checks whether the control samples are specified
        samples = ''
        samples = samples + " ".join(['-t'] + input.ChIP)

        if hasattr(input, 'Cont'):
            samples = samples + " ".join([' -c'] + input.Cont)

        # define the genome size based on sum of length of chroms
        with open(input.chrlen[0]) as f:
            chrsize = [ int(line[1]) for line in csv.reader(f,delimiter="\t")]
            gsize = sum( chrsize )

        command = " ".join(
        [params.macs2, 'callpeak',
        samples,
        '--gsize',str(gsize),
        '--outdir', params.outpath,
        '-n', params.name,
        join_params("macs2", PARAMS, params_macs),
        '2>', log.logfile
        ])
        shell(command)

# ----------------------------------------------------------------------------- #
def get_macs_output(wc):
    suffix = get_macs2_suffix(wc.name, CUSTOM_PARAMS)
    peak_name = os.path.join(PATH_PEAK, wc.name, wc.name + '_peaks.' + suffix)
    return(peak_name)

rule sort_peak:
    input:
        get_macs_output
    output:
        outfile = os.path.join(PATH_PEAK, "{name}", "{name}_qsort.bed")
    params:
        sort = SOFTWARE['sort']['executable'],
        cut  = SOFTWARE['cut']['executable']
    message:"""
            Sorting peak:
                input : {input}
                output: {output}
        """
    shell:"""
        {params.sort} -r -k9 -n {input} | {params.cut} -f1-9 > {output.outfile}
    """
