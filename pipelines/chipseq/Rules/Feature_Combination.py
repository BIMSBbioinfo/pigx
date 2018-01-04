# ----------------------------------------------------------------------------- #
rule feature_combination:
    input:
        features   = peak_files,
        annotation = rules.prepare_annotation.output.outfile,
        bw         = BW
    output:
        outfile = os.path.join(PATH_FEATURE,'Feature_Combination.tsv')
    params:
        threads     = 1,
        mem         = '8G',
        annotation  = ANNOTATION,
        outpath     = PATH_FEATURE,
        scriptdir   = SCRIPT_PATH,
        tss_up      = 2000,
        R           = SOFTWARE['R']
    log:
        log = os.path.join(PATH_LOG, 'feature_combination.log')
    message:"""
            Running: feature_combination:
                output: {output.outfile}
            """
    script:
        os.path.join(SCRIPT_PATH, 'Feature_Combination.R')

        # command = ' '.join([
        #     params.R,
        #     os.path.join(SCRIPT_PATH, 'Feature_Combination.R'),
        #     '--features',   list_to_string(input.features),
        #     '--annotation', list_to_string(params.annotation.values()),
        #     '--bw',         list_to_string(input.bw),
        #     '--scriptdir',  params.scriptdir,
        #     '--outpath',    output.outfile,
        #     '--tss.up',     str(params.tss_up)
        # ])
        # shell(command)


# ---------------------------------------------------------------------------- #
# given a list of lists, returns a flattened version
def flatten(l):
    out = []
    for item in l:
        if isinstance(item, (list, tuple)):
            out.extend(flatten(item))
        else:
            out.append(item)
    return out

def list_to_string(l):
    s = ' '.join(flatten(l))
    return s
