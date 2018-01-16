# ---------------------------------------------------------------------------- #
def get_feature_combination_infiles(name):
    infiles = []
    for file in config['feature_combination'][str(name)]:
        # This check is temporary. Once Check_Config is updated, can be removed.
        if(not file in set(PEAK_NAME_LIST.keys())):
            quit()
    
        infiles = infiles + [PEAK_NAME_LIST[file]]
    
    return(infiles)

# ----------------------------------------------------------------------------- #
rule feature_combination:
    input:
        unpack(get_feature_combination_infiles)
    output:
        outfile  = os.path.join(PATH_RDS_FEATURE,'{name}_FeatureCombination.rds')
    params:
        threads     = 1,
        mem         = '8G',
        scriptdir   = SCRIPT_PATH,
    log:
        log = os.path.join(PATH_LOG, '{name}_feature_combination.log')
    message:
        """
            Running: feature_combination:
                features:   {input}
                output:     {output.outfile}
            """
    script:
        os.path.join(SCRIPT_PATH, 'Feature_Combination.R')
