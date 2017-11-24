#----------------------------------------------------------------------------- #
rule extract_signal_annotation:
        input:
            annotation = ANNOTATION,
            wig        = rules.bam2bigWig
        output:
            outfile    = os.path.join(PATH_RDS_ANALYSIS,'{name}','{name}.Extract_Signal_Annotation.rds')
        params:
            threads     = 1,
            mem         = '16G',
            peakname    = '{name}',
            scriptdir   = SCRIPT_PATH
        log:
            log = os.path.join(PATH_LOG, 'prepare_annotation.log')
        message:"""
                Running: extract_signal_annotation:
                    wig: {input:wig}
                    output: {output.outfile}
            """
        script:
            os.path.join(SCRIPT_PATH, 'Extract_Signal_Annotation.R')




#----------------------------------------------------------------------------- #
rule extract_signal_peaks:
        input:
            bed        = rules.sort_peak.outfile,
            wig        = rules.bam2bigWig
        output:
            outfile    = os.path.join(PATH_PEAK,'{name}','{name}.Extract_Signal_Peaks.rds')
        params:
            threads     = 1,
            mem         = '16G',
            expand_peak = config['params']['rule extract_signal']['expand_peak'],
            bin_num     = config['params']['rule extract_signal']['bin_num'],
            peakname    = '{name}',
            scriptdir   = SCRIPT_PATH
        log:
            log = os.path.join(PATH_LOG, 'prepare_annotation.log')
        message:"""
                Running: extract_signal_peaks:
                    bed: {input:bed}
                    wig: {input:wig}
                    output: {output.outfile}
            """
        script:
            os.path.join(SCRIPT_PATH, 'Extract_Signal_Peaks.R')

#----------------------------------------------------------------------------- #
# Annotate Peaks
rule annotate_peaks:
        input:
            annotation = ANNOTATION,
            peaks      = rules.sort_peak.outfile,
        output:
            outfile    = os.path.join(PATH_PEAK,'{name}','{name}.Annotate_Peaks.rds')
        params:
            threads     = 1,
            mem         = '16G',
            peakname    = '{name}',
            scriptdir   = SCRIPT_PATH
        log:
            log = os.path.join(PATH_LOG, 'prepare_annotation.log')
        message:"""
                Running: annotate_peaks:
                    bed: {input:bed}
                    wig: {input:wig}
                    output: {output.outfile}
            """
        script:
            os.path.join(SCRIPT_PATH, 'Annotate_Peaks.R')
