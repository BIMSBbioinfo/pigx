# ----------------------------------------------------------------------------- #
rule parse_bowite2_log:
    input:
        bamfiles = BAMFILES_LIST
    output:
        outfile  = os.path.join(PATH_RDS, "BowtieLog.rds")
    params:
        path_log     = PATH_LOG,
        threads      = 1,
        mem          = '1G',
        scriptdir    = SCRIPT_PATH,
        Rscript      = SOFTWARE['Rscript']['executable'],
        genome_types = list(GENOME_TYPES.values())
    log:
       logfile = os.path.join(PATH_LOG, 'parse_bowite2_log.log')
    message:
        """
            Running: parse_bowite2_log:
                output:  {output.outfile}
        """
    run:
        RunRscript(input, output, params, log.logfile, 'Parse_Bowtie2log.R')
