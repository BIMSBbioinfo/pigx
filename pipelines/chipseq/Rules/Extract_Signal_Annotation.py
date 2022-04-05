#----------------------------------------------------------------------------- #
rule extract_signal_annotation:
        input:
            annotation = rules.prepare_annotation.output.outfile,
            wig        = os.path.join(PATH_MAPPED, GENOME_TYPES['Main'], "{name}", "{name}.bw")
        output:
            outfile    = os.path.join(PATH_RDS_TEMP,'{name}','{name}.Extract_Signal_Annotation.rds')
        params:
            peakname       = '{name}',
            scriptdir      = SCRIPT_PATH,
            Rscript        = SOFTWARE['Rscript']['executable'],
            number_of_bins = PARAMS['extract_signal']['number_of_bins']
        log:
            logfile = os.path.join(PATH_LOG, '{name}', '{name}.extract_signal_annotation.log')
        message:"""
                Running: extract_signal_annotation:
                    annot:  {input.annotation}
                    wig:    {input.wig}
                    output: {output.outfile}
            """
        run:
            RunRscript(input, output, params, log.logfile, 'Extract_Signal_Annotation.R')
