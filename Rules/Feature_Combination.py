# ----------------------------------------------------------------------------- #
rule feature_combination:
    input:
        features = peak_files,
        bw = BW
    output:
        outfile = os.path.join(PATH_FEATURE,'Feature_Combination.tsv')
    params:
        threads     = 1,
        mem         = '8G',
        annotation  = ANNOTATION,
        outpath     = PATH_FEATURE,
        scriptdir   = SCRIPT_PATH
    log:
        log = os.path.join(PATH_LOG, 'feature_combination.log')
    message:"""
            Running: feature_combination:
                output: {output.outfile}
            """
    script:
        os.path.join(SCRIPT_PATH, 'Feature_Combinaton.R')

