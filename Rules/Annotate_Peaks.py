#----------------------------------------------------------------------------- #
# Annotate Peaks
rule annotate_peaks:
        input:
            annotation = rules.prepare_annotation.outfile.output,
            peaks      = rules.sort_peak.outfile,
        output:
            outfile    = os.path.join(PATH_PEAK,'{name}','{name}.Annotate_Peaks.rds')
        params:
            threads     = 1,
            mem         = '16G',
            peakname    = '{name}',
            scriptdir   = SCRIPT_PATH,
            Rscript     = SOFTWARE['Rscript']['executable']
        log:
            log = os.path.join(PATH_LOG, 'prepare_annotation.log')
        message:"""
                Running: annotate_peaks:
                    annot: {input.annotation}
                    peaks: {input.peaks}
                    output: {output.outfile}
            """
        run:
            RunRscript(input, output, params, 'Annotate_Peaks.R')
