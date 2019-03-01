# ----------------------------------------------------------------------------- #
rule chipqc:
    input:
        bamfile              = os.path.join(PATH_MAPPED, '{genome_type}', "{name}", "{name}" + BAM_SUFFIX),
        index                = os.path.join(PATH_MAPPED, '{genome_type}', "{name}", "{name}" + BAM_SUFFIX + ".bai"),
        logfile              = os.path.join(PATH_RDS, "BowtieLog.rds"),
        nucleotide_frequency = os.path.join(PATH_INDEX, '{genome_type}','{genome}.NucleotideFrequency.GRanges.rds'),
        annotation           = rules.prepare_annotation.output.outfile
    output:
        outfile  = os.path.join(PATH_RDS_CHIPQC, '{genome_type}' ,"{name}_{genome}_ChIPQC.rds")
    params:
        use_longest_chr = PARAMS['chipqc']['use_longest_chr'],
        sample_name     = "{name}",
        library_type    = lambda wc: get_library_type(wc.name),
        scriptdir       = SCRIPT_PATH,
        Rscript         = SOFTWARE['Rscript']['executable'],
        threads         = config['execution']['rules']['chipqc']['threads'],
        shift_window    = SHIFT_WINDOW,
        chrM_givenName  = PARAMS['chipqc']['chrM_givenName'],
        discard_chrM    = PARAMS['chipqc']['discard_chrM']
    log:
        logfile = os.path.join(PATH_LOG, '{name}_{genome}_ChIPQC.log')
    message:
        """
            Running: ChIPQC:
                bamfile: {input.bamfile}
                output:  {output.outfile}
            """
    run:
        RunRscript(input, output, params, log.logfile, 'ChIPQC.R')
