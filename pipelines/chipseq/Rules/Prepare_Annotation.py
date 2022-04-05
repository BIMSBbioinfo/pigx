#----------------------------------------------------------------------------- #
rule link_annotation:
        input:
            annotation = ANNOTATION,
        output:
            outfile = os.path.join(PATH_ANNOTATION, 'GTF_Link.gtf')
        message:"""
                Running: link_annotation:
                    output: {output.outfile}
            """
        run:

            try: 
                os.symlink(input.annotation, output.outfile)
            except FileExistsError:
                if os.path.islink(output.outfile):
                    if not os.readlink(output.outfile) == input.annotation:
                        os.symlink(input.annotation, output.outfile)


#----------------------------------------------------------------------------- #
rule prepare_annotation:
        input:
            gtf_path = rules.link_annotation.output.outfile
        output:
            outfile = os.path.join(PATH_ANNOTATION, 'Processed_Annotation.rds')
        params:
            scriptdir    = SCRIPT_PATH,
            width_params = PARAMS['width_params'],
            Rscript      = SOFTWARE['Rscript']['executable']
        log:
            logfile = os.path.join(PATH_LOG, 'prepare_annotation.log')
        message:"""
                Running: prepare_annotation:
                    output: {output.outfile}
            """
        run:
            RunRscript(input, output, params, log.logfile, 'Prepare_Annotation.R')

