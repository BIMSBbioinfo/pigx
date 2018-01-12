#----------------------------------------------------------------------------- #
rule extract_signal_annotation:
        input:
            annotation = rules.prepare_annotation.output.outfile,
            wig        = rules.bam2bigWig.output.outfile
        output:
            outfile    = os.path.join(PATH_RDS_TEMP,'{name}','{name}.Extract_Signal_Annotation.rds')
        params:
            threads     = 1,
            mem         = '16G',
            peakname    = '{name}',
            scriptdir   = SCRIPT_PATH
        log:
            log = os.path.join(PATH_LOG, 'prepare_annotation.log')
        message:"""
                Running: extract_signal_annotation:
                    annot:  {input.annotation}
                    wig:    {input.wig}
                    output: {output.outfile}
            """
        script:
            os.path.join(SCRIPT_PATH, 'Extract_Signal_Annotation.R')




#----------------------------------------------------------------------------- #
rule extract_signal_peaks:
        input:
            bed        = rules.sort_peak.output.outfile,
            wig        = rules.bam2bigWig.output.outfile
        output:
            outfile    = os.path.join(PATH_RDS_TEMP,'{name}','{name}.Extract_Signal_Peaks.rds')
        params:
            threads     = 1,
            mem         = '16G',
            expand_peak = config['params']['extract_signal']['expand_peak'],
            bin_num     = config['params']['extract_signal']['bin_num'],
            peakname    = '{name}',
            scriptdir   = SCRIPT_PATH
        log:
            log = os.path.join(PATH_LOG, 'prepare_annotation.log')
        message:"""
                Running: extract_signal_peaks:
                    bed: {input.bed}
                    wig: {input.wig}
                    output: {output.outfile}
            """
        script:
            os.path.join(SCRIPT_PATH, 'Extract_Signal_Peaks.R')

# #----------------------------------------------------------------------------- #
# Annotate Peaks
rule annotate_peaks:
        input:
            annotation = rules.prepare_annotation.output.outfile,
            peaks      = rules.sort_peak.output.outfile
        output:
            outfile    = os.path.join(PATH_RDS_TEMP,'{name}','{name}.Annotate_Peaks.rds')
        params:
            threads     = 1,
            mem         = '16G',
            peakname    = '{name}',
            scriptdir   = SCRIPT_PATH
        log:
            log = os.path.join(PATH_LOG, 'annotate_peaks.log')
        message:"""
                Running: annotate_peaks:
                    annot: {input.annotation}
                    peaks: {input.peaks}
                    output: {output.outfile}
            """
        script:
            os.path.join(SCRIPT_PATH, 'Annotate_Peaks.R')
