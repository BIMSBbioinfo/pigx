# ---------------------------------------------------------------------------- #
def get_feature_combination_infiles(wc):
    infiles = []
    for file in config['feature_combination'][wc.name]:
        infiles = infiles + [PEAK_NAME_LIST[file]]

    infiles = dict(zip(infiles, infiles))
    return(infiles)

# ----------------------------------------------------------------------------- #
rule feature_combination:
    input:
        unpack(get_feature_combination_infiles)
    output:
        outfile  = os.path.join(PATH_RDS_FEATURE,'{name}_FeatureCombination.rds'),
        txtfile  = os.path.join(PATH_RDS_FEATURE,'{name}_FeatureCombination.txt')
    params:
        scriptdir   = SCRIPT_PATH,
        Rscript     = SOFTWARE['Rscript']['executable']
    log:
        logfile = os.path.join(PATH_LOG, '{name}_feature_combination.log')
    message:
        """
            Running: feature_combination:
                features:   {input}
                output:     {output}
            """
    run:
        RunRscript(input, output, params, log.logfile, 'Feature_Combination.R')
