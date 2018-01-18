# ----------------------------------------------------------------------------- #
rule chipqc:
    input:
        peakfiles =
        bamfiles  =
        config    = config
    output:
        outfile  = os.path.join(PATH_RDS_TEMP,'Peaks_ChIPQC.rds')
    params:
        threads     = 1,
        mem         = '8G',
        scriptdir   = SCRIPT_PATH,
    log:
        log = os.path.join(PATH_LOG, 'Peaks_ChIPQC.log')
    message:
        """
            Running: feature_combination:
                output:     {output.outfile}
            """
    run:
        RunRscript(input, output, params, 'ChIPQC.R')
