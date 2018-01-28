# ----------------------------------------------------------------------------- #
rule multiqc:
    input:
        bowtie_output   = expand(os.path.join(PATH_MAPPED, "{name}", "{name}.sorted.bam"), name=NAMES),
        fastqc_output   = FASTQC
    output:
        os.path.join(PATH_REPORTS, "multiqc.html")
    params:
        path_lot    = PATH_LOG,
        path_fastqc = PATH_QC,
        threads     = 1,
        mem         = '8G',
        # multiqc = SOFTWARE['multiqc']['executable']
        multiqc     = 'multiqc'
    log:
        os.path.join(PATH_LOG, 'multiqc.log')
    message:"""
            multiqc
        """
    shell: """
        {params.multiqc} -o {output} {params.path_log} {params.path_fastqc} >> {log} 2>&1
    """
