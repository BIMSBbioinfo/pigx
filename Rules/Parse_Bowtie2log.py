# ----------------------------------------------------------------------------- #
rule parse_bowite2_log:
    input:
        bamfiles = expand(os.path.join(PATH_MAPPED, "{name}", "{name}" + BAM_SUFFIX + ".bai"), name=NAMES)
    output:
        outfile  = os.path.join(PATH_RDS, "BowtieLog.rds")
    params:
        logfiles = expand(os.path.join(PATH_LOG, "{name}.bowtie2.log"), name = NAMES),
        threads   = 1,
        mem       = '1G',
        scriptdir = SCRIPT_PATH,
        Rscript   = SOFTWARE['Rscript']['executable']
    log:
       logfile = os.path.join(PATH_LOG, 'parse_bowite2_log.log')
    message:
        """
            Running: parse_bowite2_log:
                output:  {output.outfile}
        """
    run:
        RunRscript(input, output, params, log.logfile, 'Parse_Bowtie2log.R')
