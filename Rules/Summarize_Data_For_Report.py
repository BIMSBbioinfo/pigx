#----------------------------------------------------------------------------- #
rule summarize_data_for_report:
        input:
            REPORT_INPUT
        output:
            outfile    = SUMMARIZED_DATA_FOR_REPORT
        params:
            analysis_path   = PATH_RDS,
            # analysis names are parts of the analysis that are executed
            analysis_names  = ANALISYS_NAMES,
            # report chunks are parts of the analysis that are implemented
            report_chunks   = list(REPORT_CHUNKS.values()),
            threads         = 1,
            mem             = '32G',
            script_path     = SCRIPT_PATH,
            Rscript         = SOFTWARE['Rscript']['executable']
        log:
            logfile = os.path.join(PATH_LOG, 'Summarize_Data_For_Report.log')
        message:"""
                Running: Summarize_Data_For_Report:
                    output: {output.outfile}
            """
        run:
            RunRscript(input, output, params, log.logfile, 'Summarize_Data_For_Report.R')
