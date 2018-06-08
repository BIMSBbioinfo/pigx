#----------------------------------------------------------------------------- #
# Annotate Peaks
rule annotate_peaks:
        input:
            annotation = rules.prepare_annotation.output.outfile,
            peaks      = rules.sort_peak.output.outfile,
        output:
            outfile    = os.path.join(PATH_RDS_TEMP,'{name}','{name}.Annotate_Peaks.rds')
        params:
            peakname    = '{name}',
            scriptdir   = SCRIPT_PATH,
            Rscript     = SOFTWARE['Rscript']['executable']
        log:
            logfile = os.path.join(PATH_LOG, '{name}.Annotate_Peaks.log')
        message:"""
                Running: annotate_peaks:
                    annot: {input.annotation}
                    peaks: {input.peaks}
                    output: {output.outfile}
            """
        run:
            RunRscript(input, output, params, log.logfile, 'Annotate_Peaks.R')
