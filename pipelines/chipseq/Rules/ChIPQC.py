# ----------------------------------------------------------------------------- #
rule chipqc:
    input:
        bamfile              = os.path.join(PATH_MAPPED, "{name}", "{name}" + BAM_SUFFIX),
        index                = os.path.join(PATH_MAPPED, "{name}", "{name}" + BAM_SUFFIX + ".bai"),
        logfile              = rules.parse_bowite2_log.output.outfile,
        nucleotide_frequency = rules.extract_nucleotide_frequency.output.outfile,
        annotation           = rules.prepare_annotation.output.outfile
    output:
        outfile  = os.path.join(PATH_RDS_CHIPQC, "{name}_ChIPQC.rds")
    params:
        use_longest_chr = 'TRUE',
        sample_name     = "{name}",
        library_type    = lambda wc: get_library_type(wc.name),
        scriptdir       = SCRIPT_PATH,
        Rscript         = SOFTWARE['Rscript']['executable'],
        threads         = config['execution']['rules']['chipqc']['threads']
    log:
        logfile = os.path.join(PATH_LOG, '{name}_ChIPQC.log')
    message:
        """
            Running: ChIPQC:
                bamfile: {input.bamfile}
                output:  {output.outfile}
            """
    run:
        RunRscript(input, output, params, log.logfile, 'ChIPQC.R')
