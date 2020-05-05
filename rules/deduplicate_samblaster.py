# ==========================================================================================
# Deduplication by samblaster 
#
# https://github.com/GregoryFaust/samblaster
#
rule samblaster_markdup_sort:
    input:
        DIR_mapped+"{sample}.bwameth.bam"
    output:
        bam = DIR_sorted+"{sample}.bwameth.sorted.markdup.bam",
        index = DIR_sorted+"{sample}.bwameth.sorted.markdup.bam.bai"
    params:
        threads=config['execution']['rules']['samblaster_markdup_sort']['threads'],
        memory=config['execution']['rules']['samblaster_markdup_sort']['memory'],
        tmpdir=DIR_sorted+"{sample}/"
    log:
        DIR_sorted+"{sample}_markdups.log"
    message: fmt("Deduplicating reads with samblaster for sample {{sample}}")
    shell:
        nice("samtools", 
        ["view -h {input}"," | ", 
        tool("samblaster"),toolArgs("samblaster"),"2> {log}","|",
        tool("samtools"),"sort","-T={params.tmpdir}",
         "-o {output.bam}", "-@ {params.threads}", 
         "-m {params.memory}", "-l 9","2> {log}",";",
         tool("samtools"),"index {output.bam}"],("{log}"))
