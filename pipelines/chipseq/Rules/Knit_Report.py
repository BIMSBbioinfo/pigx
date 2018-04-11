#----------------------------------------------------------------------------- #
rule knit_report:
        input:
            infile = rules.summarize_data_for_report.output.outfile
        output:
            outfile    = REPORT
        params:
            report_template = REPORT_TEMPLATE,
            analysis_path   = PATH_RDS,
            # analysis names are parts of the analysis that are executed
            analysis_names  = ANALISYS_NAMES,
            # report chunks are parts of the analysis that are implemented
            report_chunks   = list(REPORT_CHUNKS.values()),
            # width params define the look of the profile plots
            width_params    = PARAMS['width_params'],
            logo            = os.path.join(config['locations']['pkgdatadir'], "images/Logo_PiGx.png") if os.getenv("PIGX_UNINSTALLED") else os.path.join(config['locations']['pkgdatadir'], "Logo_PiGx.png"), 
            threads         = config['execution']['rules']['knit_report']['threads'],
            mem             = config['execution']['rules']['chipqc']['memory'],
            script_path     = SCRIPT_PATH,
            Rscript         = SOFTWARE['Rscript']['executable']
        log:
            logfile = os.path.join(PATH_LOG, 'knit_report.log')
        message:"""
                Running: knit_report:
                    input:  {input.infile}
                    output: {output.outfile}
            """
        run:
            RunRscript(input, output, params, log.logfile, 'Knit_Report.R')
