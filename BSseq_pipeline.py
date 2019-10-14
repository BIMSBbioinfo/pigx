# PiGx BSseq Pipeline.
#
# Copyright © 2017, 2018 Bren Osberg <Brendan.Osberg@mdc-berlin.de>
# Copyright © 2017 Alexander Gosdschan <alexander.gosdschan@mdc-berlin.de>
# Copyright © 2019 Alexander Blume <alexander.blume@mdc-berlin.de>
# Copyright © 2017 Katarzyna Wreczycka <katwre@gmail.com>
# Copyright © 2017, 2018 Ricardo Wurmus <ricardo.wurmus@mdc-berlin.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os

#--- DEFINE OUTPUT DIRECTORIES TO BE PRODUCED 
OUTDIR = config['locations']['output-dir']                      #--- current work dir (important for rmarkdown)

DIR_scripts   = os.path.join(config['locations']['pkglibexecdir'], 'scripts/')
DIR_templates = os.path.join(OUTDIR, 'pigx_work/report_templates/')

DIR_diffmeth    =  os.path.join(OUTDIR, '09_differential_methylation/' )
DIR_seg         =  os.path.join(OUTDIR, '08_segmentation/' )
DIR_bigwig      =  os.path.join(OUTDIR, '07_bigwig_files/')
DIR_methcall    =  os.path.join(OUTDIR, '06_methyl_calls/' )
DIR_sorted      =  os.path.join(OUTDIR, '05_sorting_deduplication/' )
DIR_mapped      =  os.path.join(OUTDIR, '04_mapping/' )
DIR_posttrim_QC =  os.path.join(OUTDIR, '03_posttrimming_QC/' )
DIR_trimmed     =  os.path.join(OUTDIR, '02_trimming/' )
DIR_rawqc       =  os.path.join(OUTDIR, '01_raw_QC/' )

DIR_final       = os.path.join(OUTDIR, "Final_Reports/")


#--- DEFINE PATHS AND FILE NAMES:

PATHIN     = os.path.join(OUTDIR, "pigx_work/input/")           # location of the data files to be imported (script creates symbolic link)
GENOMEPATH = config['locations']['genome-dir']+"/"       # where the reference genome being mapped to is stored
ASSEMBLY   = config['general']['assembly'] # version of the genome being mapped to

## FIXME: request the genome file instead of the folder in the settings file
# this hack actually will only take the first file ending with either .fa or .fasta as reference genome
GENOMEFILE = os.path.join(GENOMEPATH, [file for file in os.listdir(GENOMEPATH) if file.endswith((".fasta",".fa"))][0])

#--- COMPATIBILTY VARIABLES FOR KASIAS WGBS PIPELINE:
SUBSET_READS = False
NOTRIMMING   = False
USEBWAMETH   = config['general']['use_bwameth']

# include function definitions and extra rules
include   : os.path.join(config['locations']['pkglibexecdir'], 'scripts/func_defs.py')
validate_config(config)

#--- LIST THE OUTPUT FILES TO BE PRODUCED: 

# Below is a mapping of rule names to the expected output files they
# produce.  The desired output files are specified in
# "OUTPUT_FILES".  A different set of output files can be
# selected to run fewer rules.

targets = {
    # rule to print all rule descriptions
    'help': {
        'description': "Print all rules and their descriptions.",
        'files': []
    },

    # This is an expensive one-time rule to prepare the genome.
    'genome-prep': {
        'description': "Convert reference genome into Bisulfite analogue.",
        'files': [
            GENOMEPATH+"Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
            GENOMEPATH+"Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa"
        ]
    },

    'raw-qc': {
        'description': "Perform raw quality control.",
        'files': files_for_sample(list_files_rawQC)
    },

    # This rule is always executed, as trimming is a prerequisite for
    # subsequent rules
    'trimgalore': {
        'description': "Trim the reads.",
        'files': files_for_sample(list_files_TG)
    },

    # fastQC output files are not needed downstream and need to be
    # called explicitly.
    'posttrim-qc': {
        'description': "Perform quality control after trimming.",
        'files': files_for_sample(list_files_posttrim_QC)
    },

    'mapping': {
        'description': "Align and map reads with Bismark.",
        'files': files_for_sample(list_files_bismark)
    },

    'mapping-bwameth': {
        'description': "Align and map reads with BWA-Meth.",
        'files': files_for_sample(list_files_bwameth)
    },

    'sorting': {
        'description': "Sort bam files.",
        'files': files_for_sample(list_files_sortbam)
    },

    'deduplication': {
        'description': "Deduplicate Bismark bam files.",
        'files': files_for_sample(list_files_dedupe)
    },

    'markduplicates': {
        'description': "Mark duplicates and sort BWA-Meth bam files.",
        'files': files_for_sample(list_files_markdup)
    },

    'bwameth-mapping-stats': {
        'description': "Get stats on BWA-Meth bam files.",
        'files': files_for_sample(list_files_bwamethMappingStats)
    },

    'methyl-extraction': {
        'description': "Process bwameth bam files using methylDackel.",
        'files': files_for_sample(list_files_methyldackel_extract)
    },

     # TODO: had to add this part to call bam_methCall for diff meth rule
    'methyl-calling': {
        'description': "Process bismark bam files.",
        'files': files_for_sample(bam_processing)
    },

    'maketabix-methyldackel': {
        'description': "Generate Tabix files from methylDackel files.",
        'files': files_for_sample(list_files_maketabix_methyldackel)
    },

    'export-bigwig': {
        'description': "export bigwig files from tabix files for visualization",
        'files': files_for_sample(bigwig_exporting_bwameth)
    },

    'bigwig': {
        'description': "export bigwig files to separate folder for visualization",
        'files': files_for_sample(bigwig_exporting_bismark)
    },

    'segmentation': {
        'description': "Segmentation of the methylation signal.",
        'files': files_for_sample(methSeg_bismark)
    },

    'segmentation-bwameth': {
        'description': "Segmentation of the methylation signal.",
        'files': files_for_sample(methSeg_bwameth)
    },
    'diffmeth': {
        'description': "Perform differential methylation calling.",
        'files': [ [DIR_diffmeth+"_".join(x)+
                    dedupe_tag(config["SAMPLES"][get_sampleids_from_treatment(x[0])[0]]['Protocol']) + "_diffmeth.RDS"] for x in config["general"]["differential-methylation"]["treatment-groups"] if x ]
    },

    'diffmeth-report': {
        'description': "Produce a comprehensive report for differential methylation.",
        'files':[ [DIR_final+"diffmeth-report."+"vs".join(x)+".html"] for x in config["general"]["differential-methylation"]["treatment-groups"] if x ]
    },

    'final-report': {
        'description': "Produce a comprehensive report.  This is the default target.",
        'files': files_for_sample(list_final_reports)
    }
}

TREAT_GROUPS=config["general"]["differential-methylation"]["treatment-groups"]
if TREAT_GROUPS:
  selected_targets_default = ['final-report', 'diffmeth-report', 'bigwig']
else:
  selected_targets_default = ['final-report', 'bigwig']

# Selected output files from the above set.
selected_targets = config['execution']['target'] or selected_targets_default
if USEBWAMETH: 
    selected_targets.append('mapping-bwameth')


# FIXME: the list of files must be flattened twice(!).  We should make
# sure that the targets really just return simple lists.
from itertools import chain
OUTPUT_FILES = list(chain.from_iterable(chain.from_iterable([targets[name]['files'] for name in selected_targets])))


# ==============================================================================================================
#
#                                         BEGIN RULES
#
# rules are separated by "==" bars into pairs for paired-end and single-end (subdivided by smaller "--" dividers)
# rules are (generally) presented in hierarchical order of dependency (i.e. last to first)
# ===============================================================================================================


rule all:
    input:
        OUTPUT_FILES

rule help:
    run:
        for key in sorted(targets.keys()):
            print('{}:\n  {}'.format(key, targets[key]['description']))

# Record any existing output files, so that we can detect if they have
# changed.
expected_files = {}
onstart:
    if OUTPUT_FILES:
        for name in OUTPUT_FILES:
            if os.path.exists(name):
                expected_files[name] = os.path.getmtime(name)

# Print generated target files.
onsuccess:
    if OUTPUT_FILES:
        # check if any existing files have been modified
        generated = []
        for name in OUTPUT_FILES:
            if name not in expected_files or os.path.getmtime(name) != expected_files[name]:
                generated.append(name)
        if generated:
            print("The following files have been generated:")
            for name in generated:
                print("  - {}".format(name))


# ==========================================================================================
# Generate the final report for individual samples:

rule final_report:
    input:
        rdsfile     = os.path.join(DIR_methcall,"{prefix}_methylRaw.RDS"),
        callFile    = os.path.join(DIR_methcall,"{prefix}_CpG.txt"),
        grfile      = os.path.join(DIR_seg,"{prefix}_meth_segments_gr.RDS"),
        bedfile     = os.path.join(DIR_seg,"{prefix}_meth_segments.bed"),
        template            = os.path.join(DIR_templates,"index.Rmd"),
        bigwigFile          = os.path.join(DIR_bigwig,   "{prefix}.bw"),
        chrom_seqlengths    = os.path.join(DIR_mapped,"Refgen_"+ASSEMBLY+"_chromlengths.csv")
    output:
        report        = os.path.join(DIR_final, "{prefix}_{assembly}_final.html")
    params:
        ## absolute path to bamfiles
        Samplename  = lambda wc: get_fastq_name( wc.prefix ),
        chrom_seqlengths  = os.path.join(DIR_mapped,"Refgen_"+ASSEMBLY+"_chromlengths.csv"),
        source_dir  = config['locations']['input-dir'],
        out_dir     = OUTDIR,
        bigwigFile  = os.path.join(DIR_bigwig,   "{prefix}.bw"),
        inBam       = os.path.join(OUTDIR, DIR_sorted,"{prefix}.bam"),
        assembly    = ASSEMBLY,
        mincov         = int(config['general']['methylation-calling']['minimum-coverage']),
        minqual        = int(config['general']['methylation-calling']['minimum-quality']),
        TSS_plotlength = int(config['general']['reports']['TSS_plotlength']),
        ## absolute path to output folder in working dir
        methCallRDS     = os.path.join(OUTDIR,DIR_methcall,"{prefix}_methylRaw.RDS"),
        methSegGR       = os.path.join(OUTDIR,DIR_seg,"{prefix}_meth_segments_gr.RDS"),
        methSegBed      = os.path.join(OUTDIR,DIR_seg,"{prefix}_meth_segments.bed"),
        methSegPng      = os.path.join(OUTDIR,DIR_seg,"{prefix}_meth_segments.png"),
        genome_dir  = config['locations']['genome-dir'],
        scripts_dir = DIR_scripts,
        refGenes_bedfile  = config['general']['differential-methylation']['annotation']['refGenes_bedfile'],
        webfetch    = config['general']['differential-methylation']['annotation']['webfetch']
    log:
        os.path.join(DIR_final,"{prefix}_{assembly}_final.log")
    message: fmt("Compiling final report: Output log file for details (for cluster jobs, see pigx_work/cluster_log_files).")
    run:
        generateReport(input, output, params, log, "")



# ==========================================================================================
# Generate the final report for differential methylation between pairs of treatment values:

rule diffmeth_report:
    input:
        bedfile            = lambda wc: makeDiffMethPath(DIR_diffmeth, '_diffmeth.bed',   wc),
        RDSdiffFile        = lambda wc: makeDiffMethPath(DIR_diffmeth, '_diffmeth.RDS',   wc),
        RDSdiffFile_hyper  = lambda wc: makeDiffMethPath(DIR_diffmeth, '_diffmethhyper.RDS',  wc),
        RDSdiffFile_hypo   = lambda wc: makeDiffMethPath(DIR_diffmeth, '_diffmethhypo.RDS',   wc),
        RDSdiffFile_nonsig = lambda wc: makeDiffMethPath(DIR_diffmeth, '_diffmethnonsig.RDS', wc),
        template           = os.path.join(DIR_templates,"diffmeth.Rmd"),
        chrom_seqlengths   = os.path.join(DIR_mapped,"Refgen_"+ASSEMBLY+"_chromlengths.csv")
    output:
        report        = os.path.join(DIR_final, "diffmeth-report.{treatment}.html")
    params:
        source_dir   = config['locations']['input-dir'],
        scripts_dir  = DIR_scripts,
        diffmeth_dir = DIR_diffmeth,
        genome_dir   = config['locations']['genome-dir'],
        out_dir      = OUTDIR,
        cpgIsland_bedfile = config['general']['differential-methylation']['annotation']['cpgIsland_bedfile'],
        refGenes_bedfile  = config['general']['differential-methylation']['annotation']['refGenes_bedfile'],
        chrom_seqlengths  = os.path.join(DIR_mapped,"Refgen_"+ASSEMBLY+"_chromlengths.csv"),
        assembly    = ASSEMBLY,
        treatment  = lambda wc: str(wc.treatment).replace('vs', '_'),
        qvalue     = float(config['general']['differential-methylation']['qvalue']),
        difference = float(config['general']['differential-methylation']['difference']),
        webfetch    = config['general']['differential-methylation']['annotation']['webfetch'],
        methylDiffBed           = lambda wc: makeDiffMethPath(DIR_diffmeth,  '_diffmeth.bed', wc),
        methylDiff_file         = lambda wc: makeDiffMethPath(DIR_diffmeth,  '_diffmeth.RDS', wc),
        methylDiff_hyper_file   = lambda wc: makeDiffMethPath(DIR_diffmeth,  '_diffmethhyper.RDS', wc),
        methylDiff_hypo_file    = lambda wc: makeDiffMethPath(DIR_diffmeth,  '_diffmethhypo.RDS', wc),
        methylDiff_nonsig_file  = lambda wc: makeDiffMethPath(DIR_diffmeth,  '_diffmethnonsig.RDS', wc)
    log:
        os.path.join(DIR_final,"diffmeth-report.{treatment}.log")
    message: fmt("Compiling differential methylation report " + "for treatment " + "{wildcards.treatment}")
    run:
        generateReport(input, output, params, log, "")


# ==========================================================================================
# Perform differential methylation analysis:

rule diffmeth:
    ## paths inside input and output should be relative
    input:
        inputfiles  = diffmeth_input_function
    output:
        methylDiff_file        = os.path.join(DIR_diffmeth, "{treatment}_diffmeth.RDS"),
        methylDiff_hyper_file  = os.path.join(DIR_diffmeth, "{treatment}_diffmethhyper.RDS"),
        methylDiff_hypo_file   = os.path.join(DIR_diffmeth, "{treatment}_diffmethhypo.RDS"),
        methylDiff_nonsig_file = os.path.join(DIR_diffmeth, "{treatment}_diffmethnonsig.RDS"),
        bedfile                = os.path.join(DIR_diffmeth, '{treatment}_diffmeth.bed')
    params:
        workdir     = OUTDIR,
        scripts_dir = DIR_scripts,
        inputfiles  = diffmeth_input_function,
        sampleids   = lambda wc: get_sampleids_from_treatment(wc.treatment),
        treatment   = lambda wc: [config["SAMPLES"][sampleid]['Treatment'] for sampleid in get_sampleids_from_treatment(wc.treatment)],
        assembly    = ASSEMBLY,
        qvalue      = float(config['general']['differential-methylation']['qvalue']),
        difference  = float(config['general']['differential-methylation']['difference']),
        mincov      = int(config['general']['methylation-calling']['minimum-coverage']),
        cores       = int(config['general']['differential-methylation']['cores']),
        methylDiff_file        = os.path.join(OUTDIR, DIR_diffmeth, "{treatment}_diffmeth.RDS"),
        methylDiff_hyper_file  = os.path.join(OUTDIR, DIR_diffmeth, "{treatment}_diffmethhyper.RDS"),
        methylDiff_hypo_file   = os.path.join(OUTDIR, DIR_diffmeth, "{treatment}_diffmethhypo.RDS"),
        methylDiff_nonsig_file = os.path.join(OUTDIR, DIR_diffmeth, "{treatment}_diffmethnonsig.RDS"),
        outBed      = os.path.join(OUTDIR,DIR_diffmeth,"{treatment}_diffmeth.bed")
    log:
        os.path.join(DIR_diffmeth+"{treatment}_diffmeth.log")
    message: fmt("Calculating differential methylation.")
    shell:
        nice('Rscript', ['{DIR_scripts}/methDiff.R',
                         '--inputfiles="{params.inputfiles}"',
                         '--sampleids="{params.sampleids}"',
                         '--treatment="{params.treatment}"',
                         '--assembly={params.assembly}',
                         '--qvalue={params.qvalue}',
                         '--difference={params.difference}',
                         '--mincov={params.mincov}',
                         '--cores={params.cores}',
                         '--methylDiff_file={params.methylDiff_file}',
                         '--methylDiff_hyper_file={params.methylDiff_hyper_file}',
                         '--methylDiff_hypo_file={params.methylDiff_hypo_file}',
                         '--methylDiff_nonsig_file={params.methylDiff_nonsig_file}',
                         '--outBed={params.outBed}',
                         '--logFile={log}'])



# ==========================================================================================
# Perform segmentation on the methylome:

rule methseg:
    ## paths inside input and output should be relative
    input:
        tabixfile   =  os.path.join(DIR_methcall,"{tool}","tabix_{context}","{prefix}_{context}.txt.bgz")
    output:
        bedfile     = os.path.join(DIR_seg,"{prefix}_{context}_{tool}.meth_segments.bed")
    params:
        methSegPng  = os.path.join(DIR_seg,"{prefix}_{context}_{tool}.meth_segments.png"),
        sample_id   = "{prefix}",
        assembly    = ASSEMBLY
    log:
        os.path.join(DIR_seg,"{prefix}_{context}_{tool}.meth_segments.log")
    message: fmt("Segmenting methylation profile for {input.tabixfile}.")
    shell:
        nice('Rscript', ["{DIR_scripts}/methSeg.R",
                         "--tabix={input.tabixfile}",
                         "--outBed={output.bedfile}",
                         "--png={params.methSegPng}",
                         "--sample.id={params.sample_id}",
                         "--assembly={params.assembly}",
                         "--logFile={log}"],"{log}")


# ==========================================================================================
# Process bam files into methyl-called formats:

rule bam_methCall:
    input:
        bamfile         = os.path.join(DIR_sorted,"{prefix}.bam")
    output:
        tabixfile       = os.path.join(DIR_methcall,"methylKit","tabix_{context}","{prefix}_{context}.txt.bgz"),
        tabixfileindex  = os.path.join(DIR_methcall,"methylKit","tabix_{context}","{prefix}_{context}.txt.bgz.tbi")
    params:
        assembly        = ASSEMBLY,
        mincov          = int(config['general']['methylation-calling']['minimum-coverage']),
        minqual         = int(config['general']['methylation-calling']['minimum-quality']),
        context         = "{context}"
    log:
        os.path.join(DIR_methcall,"{prefix}_{context}_meth_calls.log")
    message: fmt("Extract methylation calls from bam file.")
    shell:
        nice('Rscript', ["{DIR_scripts}/methCall.R",
                         "--inBam={input.bamfile}",
                         "--assembly={params.assembly}",
                         "--mincov={params.mincov}",
                         "--minqual={params.minqual}",
                         "--context={params.context}",
                         "--tabix={output.tabixfile}",
                         "--logFile={log}"],"{log}")


# ==========================================================================================
# Deduplicate aligned reads from the bam file:

rule deduplication_se:
    input:
        DIR_sorted+"{sample}_se_bt2.sorted.bam"
    output:
        DIR_sorted+"{sample}_se_bt2.sorted.deduped.bam"
    params:
        bam="--bam ",
        sampath="--samtools_path " + tool('samtools')
    log:
        DIR_sorted+"{sample}_deduplication.log"
    message: fmt("Deduplicating single-end aligned reads from {input}")
    shell:
        nice('samtools', [" markdup -rs ", "{input}", "{output}"], "{log}")

#-----------------------
rule deduplication_pe:
    input:
        DIR_sorted+"{sample}_1_val_1_bt2.sorted.bam"
    output:
        DIR_sorted+"{sample}_1_val_1_bt2.sorted.deduped.bam"
    log:
        DIR_sorted+"{sample}_deduplication.log"
    message: fmt("Deduplicating paired-end aligned reads from {input}")
    shell:
        nice('samtools', [" markdup -r ", "{input}", "{output}"], "{log}")


# ==========================================================================================
# Sort the bam file by position (and carry out mate-flagging in paired-end case):

rule sortbam_se:
    input:
        DIR_mapped+"{sample}_trimmed_bismark_bt2.bam"
    output:
        DIR_sorted+"{sample}_se_bt2.sorted.bam"
    message: fmt("Sorting bam file {input}")
    shell:
        nice('samtools', ["sort", "{input}", "-o {output}"])

#-----------------------
rule sortbam_pe:
    input:
        DIR_mapped+"{sample}_1_val_1_bismark_bt2_pe.bam"
    output:
        DIR_sorted+"{sample}_1_val_1_bt2.sorted.bam"
    message: fmt("Sorting bam file {input}")
    shell:
        nice('samtools', ["sort -n ", " {input} ", " | ", tool('samtools'), " fixmate -m  - - ", " | ", tool('samtools'), " sort -o {output} " ])


# ==========================================================================================
# Align and map reads to the reference genome using bismark:

bismark_cores = str(config['tools']['bismark']['cores'])

rule bismark_align_and_map_se:
    input:
        refconvert_CT = GENOMEPATH+"Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
	refconvert_GA = GENOMEPATH+"Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa",
        fqfile = DIR_trimmed+"{sample}_trimmed.fq.gz",
        qc     = DIR_posttrim_QC+"{sample}_trimmed_fastqc.html"
    output:
        DIR_mapped+"{sample}_trimmed_bismark_bt2.bam",
        DIR_mapped+"{sample}_trimmed_bismark_bt2_SE_report.txt"
    params:
        bismark_args = config['tools']['bismark']['args'],
        genomeFolder = "--genome_folder " + GENOMEPATH,
        outdir = "--output_dir  "+DIR_mapped,
        nucCov = "--nucleotide_coverage",
        pathToBowtie = "--path_to_bowtie "+ os.path.dirname(tool('bowtie2')),
        useBowtie2   = "--bowtie2 ",
        samtools     = "--samtools_path "+ os.path.dirname(tool('samtools')),
        tempdir      = "--temp_dir " + DIR_mapped,
        cores = "--multicore " + bismark_cores
    log:
        DIR_mapped+"{sample}_bismark_se_mapping.log"
    message: fmt("Mapping single-end reads to genome {ASSEMBLY}")
    shell:
        nice('bismark', ["{params}", "{input.fqfile}"], "{log}")

rule bismark_align_and_map_pe:
    input:
        refconvert_CT = GENOMEPATH+"Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
	refconvert_GA = GENOMEPATH+"Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa",
        fin1 = DIR_trimmed+"{sample}_1_val_1.fq.gz",
        fin2 = DIR_trimmed+"{sample}_2_val_2.fq.gz",
        qc   = [ DIR_posttrim_QC+"{sample}_1_val_1_fastqc.html",
                 DIR_posttrim_QC+"{sample}_2_val_2_fastqc.html"]
    output:
        DIR_mapped+"{sample}_1_val_1_bismark_bt2_pe.bam",
        DIR_mapped+"{sample}_1_val_1_bismark_bt2_PE_report.txt"
    params:
        bismark_args = config['tools']['bismark']['args'],
        genomeFolder = "--genome_folder " + GENOMEPATH,
        outdir       = "--output_dir  "+DIR_mapped,
        nucCov       = "--nucleotide_coverage",
        pathToBowtie = "--path_to_bowtie "+ os.path.dirname(tool('bowtie2')),
        useBowtie2   = "--bowtie2 ",
        samtools     = "--samtools_path "+ os.path.dirname(tool('samtools')),
        tempdir      = "--temp_dir "+DIR_mapped,
        cores        = "--multicore "+bismark_cores
    log:
        DIR_mapped+"{sample}_bismark_pe_mapping.log"
    message: fmt("Mapping paired-end reads to genome {ASSEMBLY}.")
    shell:
        nice('bismark', ["{params}", "-1 {input.fin1}", "-2 {input.fin2}"], "{log}")


# ==========================================================================================
# Align and map reads to the reference genome using bwa-meth:

include: './rules/Align_bwameth_rules.py'

# ==========================================================================================
# Mark duplicate reads from bwa-meth alignment using picard-markduplicates like algo:

include: './rules/deduplicate_samblaster.py'

# ==========================================================================================
# extract mapping statistics like duplicate numbers and flagstats using samtools:

include: './rules/mapping_stats.py'

# ==========================================================================================
# Extract methylation counts with methylDackel, 
# create methylation bias and export to tabix 

include: './rules/preprocessing_methyldackel.py'

# ==========================================================================================
# generate bigwig from tabix

include: './rules/export_tabix_bigwig.py'

# ==========================================================================================
# Generate methyl-converted version of the reference genome, if necessary:

rule bismark_genome_preparation:
    input:
        ancient(GENOMEPATH)
    output:
        GENOMEPATH+"Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
        GENOMEPATH+"Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa"
    params:
        bismark_genome_preparation_args = config['tools']['bismark-genome-preparation']['args'],
        pathToBowtie = "--path_to_bowtie "+ os.path.dirname(tool('bowtie2')),
        useBowtie2   = "--bowtie2 ",
        verbose      = "--verbose "
    log:
        'bismark_genome_preparation_'+ASSEMBLY+'.log'
    message: fmt("Converting {ASSEMBLY} Genome into Bisulfite analogue")
    shell:
        nice('bismark-genome-preparation', ["{params}", "{input}"], "{log}")



# ==========================================================================================
# Create a csv file tabulating the lengths of the chromosomes in the reference genome:

rule tabulate_seqlengths:
    input:
        GENOMEFILE
    output:
        index       = GENOMEFILE+".fai",
        chrom_seqlengths  = os.path.join(DIR_mapped,"Refgen_"+ASSEMBLY+"_chromlengths.csv")
    message: fmt("Tabulating chromosome lengths in genome: {ASSEMBLY} for later reference.")
    shell:
        nice('samtools', 
        ['faidx {input}',";",
        tool('cut'),"-f1,2","{output.index}","> {output.chrom_seqlengths}"])


# ==========================================================================================
# Carry out post-trimming quality control

rule fastqc_after_trimming_se:
    input:
        DIR_trimmed+"{sample}_trimmed.fq.gz",
    output:
    	DIR_posttrim_QC+"{sample}_trimmed_fastqc.html",
    	DIR_posttrim_QC+"{sample}_trimmed_fastqc.zip"
    params:
        fastqc_args = config['tools']['fastqc']['args'],
        outdir = "--outdir "+DIR_posttrim_QC
    log:
   	    DIR_posttrim_QC+"{sample}_trimmed_fastqc.log"
    message: fmt("Quality checking trimmmed single-end data from {input}")
    shell:
        nice('fastqc', ["{params}", "{input}"], "{log}")

rule fastqc_after_trimming_pe:
    input:
        DIR_trimmed+"{sample}_1_val_1.fq.gz",
        DIR_trimmed+"{sample}_2_val_2.fq.gz"
    output:
    	DIR_posttrim_QC+"{sample}_1_val_1_fastqc.html",
    	DIR_posttrim_QC+"{sample}_1_val_1_fastqc.zip",
    	DIR_posttrim_QC+"{sample}_2_val_2_fastqc.zip",
        DIR_posttrim_QC+"{sample}_2_val_2_fastqc.html"
    params:
        fastqc_args = config['tools']['fastqc']['args'],
        outdir = "--outdir "+DIR_posttrim_QC
    log:
   	    DIR_posttrim_QC+"{sample}_trimmed_fastqc.log"
    message: fmt("Quality checking trimmmed paired-end data from {input}")
    shell:
        nice('fastqc', ["{params}", "{input}"], "{log}")


# ==========================================================================================
# Trim the reads for adapter-ends and quality

rule trim_reads_se:
    input:
       qc   = DIR_rawqc+"{sample}_fastqc.html",
       file = PATHIN+"{sample}.fq.gz"
    output:
       DIR_trimmed+"{sample}_trimmed.fq.gz" #---- this ALWAYS outputs .fq.qz format.
    params:
       extra      = config['tools']['trim-galore']['args'],
       outdir     = "--output_dir "+DIR_trimmed,
       phred      = "--phred33",
       gz         = "--gzip",
       cutadapt   = "--path_to_cutadapt " + tool('cutadapt'),
    log:
       DIR_trimmed+"{sample}.trimgalore.log"
    message: fmt("Trimming raw single-end read data from {input}")
    shell:
       nice('trim-galore', ["{params}", "{input.file}"], "{log}")

rule trim_reads_pe:
    input:
        qc    = [ DIR_rawqc+"{sample}_1_fastqc.html",
                  DIR_rawqc+"{sample}_2_fastqc.html"],
        files = [ PATHIN+"{sample}_1.fq.gz",
                  PATHIN+"{sample}_2.fq.gz"]
    output:
        DIR_trimmed+"{sample}_1_val_1.fq.gz", #---- this ALWAYS outputs .fq.qz format.
        DIR_trimmed+"{sample}_2_val_2.fq.gz",
    params:
        extra          = config['tools']['trim-galore']['args'],
        outdir         = "--output_dir "+DIR_trimmed,
        phred          = "--phred33",
        gz             = "--gzip",
        cutadapt       = "--path_to_cutadapt " + tool('cutadapt'),
        paired         = "--paired"
    log:
        DIR_trimmed+"{sample}.trimgalore.log"
    message:
        fmt("Trimming raw paired-end read data from {input}")
    shell:
        nice('trim-galore', ["{params}", "{input.files}"], "{log}")


# ==========================================================================================
# Perform quality control on raw data

rule fastqc_raw: #----only need one: covers BOTH pe and se cases.
    input:
        PATHIN+"{sample}.fq.gz"
    output:
        DIR_rawqc+"{sample}_fastqc.html",
        DIR_rawqc+"{sample}_fastqc.zip"
    params:
        fastqc_args = config['tools']['fastqc']['args'],
        outdir      = "--outdir "+ DIR_rawqc     # usually pass params as strings instead of wildcards.
    log:
        DIR_rawqc+"{sample}_fastqc.log"
    message: fmt("Quality checking raw read data from {input}")
    shell:
        nice('fastqc', ["{params}", "{input}"], "{log}")


