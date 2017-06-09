# ----------------------------------------------------------------------------- #
rule bowtie:
	input:
		"Fastq/{name}.fastq.gz"
	output:
		os.path.join(MAPPED_DIR, "{name}/{name}.bam")
	params:
		genome=GENOME,
		threads=8,
		mem='16G'
	log:
		os.path.join(MAPPED_DIR, "{name}/{name}.bowtie.log")
	message: """
		Mapping ...
	"""
	shell: 	"""
		zcat {input} | bowtie -p {params.threads} -S -k 1 -m 1 --best --strata {params.genome} -  2> {log} | samtools view -bhS > {output}
	"""

# ----------------------------------------------------------------------------- #
rule samtools_sort:
	input:
		os.path.join(MAPPED_DIR, "{name}/{name}.bam")
	output:
		os.path.join(MAPPED_DIR, "{name}/{name}.sorted.bam")
	params:
		threads = 4,
		mem = '16G'
	message: """
		Sorting ...
	"""
	shell: """
		samtools sort --threads {params.threads} -o {output} {input}
	"""

# ----------------------------------------------------------------------------- #
rule samtools_index:
	input:
		os.path.join(MAPPED_DIR, "{name}/{name}.sorted.bam")
	output:
		os.path.join(MAPPED_DIR, "{name}/{name}.sorted.bam.bai")
	params:
		threads = 1,
		mem = '8G'
	message: """
		Index ...
	"""
	shell:
		"samtools index {input}"


# ----------------------------------------------------------------------------- #
rule fastqc:
	input:
	"Fastq/{name}.fastq.gz"
	output:
		os.path.join(QC_PATH, "{name}.fastqc.zip")
	params:
		threads = 1,
		mem = '8G'
	message:
		"""FastQC"""
	shell: """
		fastqc --outdir  {input} --extract  -f fastq {output}
	"""

# ----------------------------------------------------------------------------- #
rule chrlen:
	input:
		os.path.join(MAPPED_DIR, "{name}/{name}.sorted.bam")
	output:
		os.path.join(MAPPED_DIR, "{name}/chrlen.txt")
	params:
		threads = 1,
		mem = '4G'
	message: """
		chrlen ...
	"""
	shell:"""
		samtools view -H {input} | grep @SQ | perl -pe 's/^@.+?://;s/LN://' > {output}
	"""

rule bam2bed_extend:
	input:
		file = os.path.join(MAPPED_DIR, "{name}/{name}.sorted.bam"),
		chrlen = rules.chrlen.output
	params:
		check = CHECK,
		threads = 1,
		mem = '16G'
	output:
		os.path.join(MAPPED_DIR, "{name}/{name}.bed")
	log:
		os.path.join(MAPPED_DIR, "{name}/{name}.bed.log")
	message: """
		bam2bed extend...
	"""
	shell: """
	bamToBed -i {input.file} > {output}
	#perl -lane "if(@F[5] eq '+'){{print join 	(\"\t\",@F[0..1],@F[1]+200,@F[3..5]);}}else{{print join (\"\t\",@F[0],@F[2]-200,@F[2],@F[3..5]);}}" > {output}
	#| {params.check} {input.chrlen} > {output}
	"""


rule bam2bigWig:
	input:
		chrlen = rules.chrlen.output,
		bedfile = rules.bam2bed_extend.output
	output:
		os.path.join(MAPPED_DIR, "{name}/{name}.bw")
	params:
		threads = 1,
		mem = '16G',
		bedgraph = os.path.join(MAPPED_DIR, "{name}/{name}.bedGraph")
	message: """
		bam2bedgraph ...
	"""
	shell: """
	genomeCoverageBed -i {input.bedfile} -g {input.chrlen} -bg > {params.bedgraph}
	wigToBigWig {params.bedgraph} {input.chrlen} {output}
	"""
# ----------------------------------------------------------------------------- #
rule makelinks:
	input:
		file = lambda wildcards: os.path.join(os.getcwd(), MAPPED_DIR, ex_nams[wildcards.ex_name], ex_nams[wildcards.ex_name] + '.bw')
	output:
		os.path.join(MAPPED_DIR, "Tracks",  "{ex_name}.bw")
	shell: """
		ln -s {input.file} {output}
	"""
