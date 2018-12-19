# ----------------------------------------------------------------------------- #
rule trim_galore_pe:
    input: 
        infile = lambda wildcards: get_fastq_input(wildcards.sample,wildcards.prefix)
    output:
        r1 = os.path.join(PATH_TRIMMED, "{sample}", "{prefix}_R1.fastq.gz"),
        r2 = os.path.join(PATH_TRIMMED, "{sample}", "{prefix}_R2.fastq.gz")
    params:
        tmp1 = lambda wildcards: os.path.join(PATH_TRIMMED, wildcards.sample, replace_fastq_ext(os.path.basename(get_fastq_input(wildcards.sample,wildcards.prefix)[0]),'_val_1.fq.gz')),
        tmp2 = lambda wildcards: os.path.join(PATH_TRIMMED, wildcards.sample, replace_fastq_ext(os.path.basename(get_fastq_input(wildcards.sample,wildcards.prefix)[1]),'_val_2.fq.gz')),
        trim_galore = SOFTWARE['trim_galore']['executable'],
        mv = SOFTWARE['mv']['executable']
    message:"""
        Trimming Fastq files:
            input : {input}
            output: {output}
    """
    log: 
        os.path.join(PATH_LOG, "{sample}",'trim_galore_{prefix}.log')
    shell:
        "{params.trim_galore} -o {PATH_TRIMMED}/{wildcards.sample} --paired {input.infile[0]} {input.infile[1]} >> {log} 2>&1 ; {params.mv} {params.tmp1} {output.r1} && {params.mv} {params.tmp2} {output.r2}"

rule trim_galore_se:
    input: 
        infile = lambda wildcards: get_fastq_input(wildcards.sample,wildcards.prefix)
    output: 
        read = os.path.join(PATH_TRIMMED,"{sample}", "{prefix}_R.fastq.gz")
    params: 
        tmp = lambda wildcards: os.path.join(PATH_TRIMMED, wildcards.sample, replace_fastq_ext(os.path.basename(get_fastq_input(wildcards.sample,wildcards.prefix)[0]),'_trimmed.fq.gz')),
        trim_galore = SOFTWARE['trim_galore']['executable'],
        mv = SOFTWARE['mv']['executable']
    message:"""
        Trimming Fastq files:
            input : {input}
            output: {output}
    """
    log: 
        os.path.join(PATH_LOG, "{sample}",'trim_galore_{prefix}.log')
    shell:
        "{params.trim_galore} -o {PATH_TRIMMED}/{wildcards.sample} {input.infile[0]} >> {log} 2>&1 ; {params.mv} {params.tmp} {output.read}"
