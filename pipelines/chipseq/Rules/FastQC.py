# ----------------------------------------------------------------------------- #
def get_fastqc_input(wc):
    return FASTQC_DICT[wc.name]['fastq']

rule fastqc:
    input:
        infile = get_fastqc_input
    output:
        outfile = os.path.join(PATH_QC, "{folder}", "{name}_fastqc.zip")
    params:
        outpath = os.path.join(PATH_QC, "{folder}"),
        fastqc  = SOFTWARE['fastqc']['executable']
    log:
        os.path.join(PATH_LOG, "{name}.fastqc.log")
    message:"""
            FastQC:
                input: {input.infile}
                output: {output.outfile}
        """
    shell: """
        {params.fastqc} --outdir {params.outpath} --extract -f fastq -t 1 {input.infile} >> {log} 2>&1;

    """
