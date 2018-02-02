#----------------------------------------------------------------------------- #
rule extract_signal_annotation:
        input:
            annotation = rules.prepare_annotation.output.outfile,
            wig        = rules.bam2bigWig.output.outfile
        output:
            outfile    = os.path.join(PATH_RDS_TEMP,'{name}','{name}.Extract_Signal_Annotation.rds')
        params:
            peakname    = '{name}',
            scriptdir   = SCRIPT_PATH,
            Rscript     = SOFTWARE['Rscript']['executable']
        log:
            log = os.path.join(PATH_LOG, 'prepare_annotation.log')
        message:"""
                Running: extract_signal_annotation:
                    annot:  {input.annotation}
                    wig:    {input.wig}
                    output: {output.outfile}
            """
        run:
            RunRscript(input, output, params, 'Extract_Signal_Annotation.R')





#----------------------------------------------------------------------------- #
rule extract_signal_peaks:
        input:
            bed        = rules.sort_peak.output.outfile,
            wig        = rules.bam2bigWig.output.outfile
        output:
            outfile    = os.path.join(PATH_RDS_TEMP,'{name}','{name}.Extract_Signal_Peaks.rds')
        params:
            expand_peak = PARAMS['extract_signal']['expand_peak'],
            bin_num     = PARAMS['extract_signal']['bin_num'],
            peakname    = '{name}',
            scriptdir   = SCRIPT_PATH,
            Rscript     = SOFTWARE['Rscript']['executable']
        log:
            log = os.path.join(PATH_LOG, 'prepare_annotation.log')
        message:"""
                Running: extract_signal_peaks:
                    bed: {input.bed}
                    wig: {input.wig}
                    output: {output.outfile}
            """
        run:
            RunRscript(input, output, params, 'Extract_Signal_Peaks.R')


# #----------------------------------------------------------------------------- #
# Annotate Peaks
rule annotate_peaks:
        input:
            annotation = rules.prepare_annotation.output.outfile,
            peaks      = rules.sort_peak.output.outfile
        output:
            outfile    = os.path.join(PATH_RDS_TEMP,'{name}','{name}.Annotate_Peaks.rds')
        params:
            peakname    = '{name}',
            scriptdir   = SCRIPT_PATH,
            Rscript     = SOFTWARE['Rscript']['executable']
        log:
            log = os.path.join(PATH_LOG, 'annotate_peaks.log')
        message:"""
                Running: annotate_peaks:
                    annot: {input.annotation}
                    peaks: {input.peaks}
                    output: {output.outfile}
            """
        run:
            RunRscript(input, output, params, 'Annotate_Peaks.R')
