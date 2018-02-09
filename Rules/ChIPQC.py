# ----------------------------------------------------------------------------- #
rule chipqc:
    input:
        bamfile              = os.path.join(PATH_MAPPED, "{name}", "{name}.sorted.bam"),
        index                = os.path.join(PATH_MAPPED, "{name}", "{name}.sorted.bam.bai"),
        logfile              = rules.parse_bowite2_log.output.outfile,
        nucleotide_frequency = rules.extract_nucleotide_frequency.output.outfile
    output:
        outfile  = os.path.join(PATH_RDS_CHIPQC, "{name}_ChIPQC.rds")
    params:
        threads         = 1,
        mem             = '8G',
        use_longest_chr = 'TRUE',
        sample_name     = "{name}",
        sample_sheet    = SAMPLE_SHEET,
        scriptdir       = SCRIPT_PATH,
        Rscript         = SOFTWARE['Rscript']['executable']
    log:
        log = os.path.join(PATH_LOG, '{name}_ChIPQC.log')
    message:
        """
            Running: ChIPQC:
                bamfile: {input.bamfile}
                output:  {output.outfile}
            """
    run:
        RunRscript(input, output, params, 'ChIPQC.R')
