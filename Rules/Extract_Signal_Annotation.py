#----------------------------------------------------------------------------- #
rule extract_signal_annotation:
        input:
            annotation = rules.prepare_annotation.output.outfile,
            wig        = rules.bam2bigWig.output.outfile
        output:
            outfile    = os.path.join(PATH_RDS_TEMP,'{name}','{name}.Extract_Signal_Annotation.rds')
        params:
            peakname       = '{name}',
            number_of_bins = 40,
            scriptdir      = SCRIPT_PATH,
            Rscript        = SOFTWARE['Rscript']['executable']
        log:
            log = os.path.join(PATH_LOG, '{name}.extract_signal_annotation.log')
        message:"""
                Running: extract_signal_annotation:
                    annot:  {input.annotation}
                    wig:    {input.wig}
                    output: {output.outfile}
            """
        run:
            RunRscript(input, output, params, 'Extract_Signal_Annotation.R')
