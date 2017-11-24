#----------------------------------------------------------------------------- #
rule prepare_annotation:
        input:
            annotation = ANNOTATION,
        output:
            outfile = os.path.join(PATH_RDS_ANALYSIS,'Processed_Annotation.rds')
        params:
            threads   = 1,
            mem       = '16G',
            scriptdir = SCRIPT_PATH
        log:
            log = os.path.join(PATH_LOG, 'prepare_annotation.log')
        message:"""
                Running: prepare_annotation:
                    output: {output.outfile}
            """
        script:
            os.path.join(SCRIPT_PATH, 'Prepare_Annotation.R')
