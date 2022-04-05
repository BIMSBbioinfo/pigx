# ---------------------------------------------------------------------------- #
rule multiqc:
    input:
        bowtie_output   = BOWTIE2,
        fastqc_output   = FASTQC,
        bamstat_output  = BAMSTATS
    output:
        os.path.join(PATH_REPORTS, "multiqc.html")
    params:
        path_log    = PATH_LOG,
        path_fastqc = PATH_QC,
        path_mapped = PATH_MAPPED,
        threads     = 1,
        mem         = '8G',
        multiqc = SOFTWARE['multiqc']['executable']
    log:
        os.path.join(PATH_LOG, 'multiqc.log')
    message:"""
            multiqc
        """
    shell: """
        {params.multiqc} -f -n {output} \
                {params.path_log} \
                {params.path_fastqc} \
                {params.path_mapped} \
                >> {log} 2>&1
    """
