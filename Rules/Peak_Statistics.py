# ----------------------------------------------------------------------------- #
rule peak_statistics:
    input:
        bamfile    = BOWTIE2,
        peaks      = QSORT,
        bwfiles    = BW,
        readnumber = rules.parse_bowite2_log.output.outfile
    output:
        outfile    = os.path.join(PATH_RDS, "Peak_Statistics.rds")
    params:
        threads         = 16,
        mem             = '16G',
        config          = SAMPLE_SHEET,
        path_mapped     = PATH_MAPPED,
        path_peak       = PATH_PEAK,
        scriptdir       = SCRIPT_PATH,
        Rscript         = SOFTWARE['Rscript']['executable'],
        peaks_resize    = 500
    log:
        log = os.path.join(PATH_LOG, 'Peak_Statistics.log')
    message:
        """
            Running: peak_statistics:
                output:  {output.outfile}
            """
    run:
        RunRscript(input, output, params, 'Peak_Statistics.R')
