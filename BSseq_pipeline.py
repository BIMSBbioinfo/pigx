# PiGx BSseq Pipeline.
#
# Copyright © 2017, 2018 Bren Osberg <Brendan.Osberg@mdc-berlin.de>
# Copyright © 2017, 2018, 2019, 2020 Alexander Blume <alexander.blume@mdc-berlin.de>
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

# include function definitions and extra rules
include   : os.path.join(config['locations']['pkglibexecdir'], 'scripts/func_defs.py')
validate_config(config)

#--- DEFINE OUTPUT DIRECTORIES TO BE PRODUCED 
OUTDIR = config['locations']['output-dir']                      #--- current work dir (important for rmarkdown)

DIR_scripts   = os.path.join(config['locations']['pkglibexecdir'], 'scripts/')
DIR_templates = os.path.join(config['locations']['pkgdatadir'], 'report_templates/')

DIR_diffmeth    =  os.path.join(OUTDIR, '09_differential_methylation/' )
DIR_seg         =  os.path.join(OUTDIR, '08_segmentation/' )
DIR_bigwig      =  os.path.join(OUTDIR, '07_bigwig_files/')
DIR_methcall    =  os.path.join(OUTDIR, '06_methyl_calls/' )
DIR_sorted      =  os.path.join(OUTDIR, '05_sorting_deduplication/' )
DIR_mapped      =  os.path.join(OUTDIR, '04_mapping/' )
DIR_posttrim_QC =  os.path.join(OUTDIR, '03_posttrimming_QC/' )
DIR_trimmed     =  os.path.join(OUTDIR, '02_trimming/' )
DIR_rawqc       =  os.path.join(OUTDIR, '01_raw_QC/' )

DIR_final       = os.path.join(OUTDIR, "Reports/")


#--- DEFINE PATHS AND FILE NAMES:


if os.getenv("PIGX_UNINSTALLED"):
    LOGOPATH = os.path.join(config['locations']['pkgdatadir'], "images/Logo_PiGx.png") 
else:
    LOGOPATH = os.path.join(config['locations']['pkgdatadir'], "Logo_PiGx.png")

BIBTEXPATH  = os.path.join(config['locations']['pkgdatadir'], "report_templates/reports.bib")

PATHIN     = os.path.join(OUTDIR, "pigx_work/input/")           # location of the data files to be imported (script creates symbolic link)
GENOMEFILE = config['locations']['genome-fasta']    # where the reference genome being mapped to is stored
GENOMEPATH = os.path.dirname(GENOMEFILE) + "/"
ASSEMBLY   = config['general']['assembly'] # version of the genome being mapped to

# FIXME: lookup and fetch now done in diffmethreport, but should they get their own rule ??
## should we do fetching at all? would be maybe more stable if we require people to download themselves.
WEBFETCH = TrueOrFalse(config['general']['differential-methylation']['annotation']['webfetch'])
CPGISLAND_BEDFILE = config['general']['differential-methylation']['annotation']['cpgIsland-bedfile']
REFGENES_BEDFILE  = config['general']['differential-methylation']['annotation']['refGenes-bedfile']

if CPGISLAND_BEDFILE and os.path.isfile(CPGISLAND_BEDFILE):
  # make path absolute
  CPGISLAND_BEDFILE = os.path.abspath(CPGISLAND_BEDFILE)
  
elif WEBFETCH:
    CPGISLAND_BEDFILE = os.path.join(OUTDIR, 'pigx_work','refGenome',"cpgIslandExt."+ASSEMBLY+".bed.gz")
    print("WARNING: Parameter 'general::differential-methylation::annotation::cpgIsland-bedfile' was not set to a valid file.\n",
    "Updating to "+CPGISLAND_BEDFILE+" since webfetch was set.\n")
    

if REFGENES_BEDFILE and os.path.isfile(REFGENES_BEDFILE):
  # make path absolute
  REFGENES_BEDFILE = os.path.abspath(REFGENES_BEDFILE)
  
elif WEBFETCH:
  REFGENES_BEDFILE = os.path.join(OUTDIR, 'pigx_work','refGenome',"knownGene."+ASSEMBLY+".bed.gz")
  print("WARNING: Parameter 'general::differential-methylation::annotation::refGenes-bedfile' was not set to a valid file.\n",
  "Updating to "+REFGENES_BEDFILE+" since webfetch was set.\n")
  
  

#--- CHOOSE PIPELINE BRANCH
USEBWAMETH = TrueOrFalse(config['general']['use_bwameth'])
USEBISMARK = TrueOrFalse(config['general']['use_bismark'])


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
    'trimming': {
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
        'description': "Process Bismark bam files.",
        'files': files_for_sample(bam_processing)
    },

    'maketabix-methyldackel': {
        'description': "Generate Tabix files from methylDackel files.",
        'files': files_for_sample(list_files_maketabix_methyldackel)
    },

    'bigwig-bwameth': {
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
    
    'unite-bwameth': {
        'description': "Unite samples for differential methylation calling.",
        'files': files_for_treatment(list_files_unite_bwameth)
    },

    'unite': {
        'description': "Unite samples for differential methylation calling.",
        'files': files_for_treatment(list_files_unite_bismark)
    },
    
    'diffmeth-bwameth': {
        'description': "Perform differential methylation calling.",
        'files': files_for_treatment(list_files_diffmeth_bwameth)
    },

    'diffmeth': {
        'description': "Perform differential methylation calling.",
        'files': files_for_treatment(list_files_diffmeth_bismark)
    },

    'diffmeth-report-bwameth': {
        'description': "Produce a comprehensive report for differential methylation.",
        'files': files_for_treatment(list_files_diffmeth_report_bwameth)
    },
    
    'diffmeth-report': {
        'description': "Produce a comprehensive report for differential methylation.",
        'files': files_for_treatment(list_files_diffmeth_report_bismark)
    },

    'final-report-bwameth': {
        'description': "Produce a comprehensive report per Sample.",
        'files': files_for_sample(list_final_reports_bwameth)
    },
    'final-report': {
        'description': "Produce a comprehensive report per Sample.",
        'files': files_for_sample(list_final_reports_bismark)
    },
    'multiqc': {
        'description': "Produce a summarized qc report for bismark branch.",
        'files': [[os.path.join(DIR_final,"multiqc","bismark_multiqc_report.html")]]
    },
    'multiqc-bwameth': {
        'description': "Produce a comprehensive report for bwameth branch.",
        'files': [[os.path.join(DIR_final,"multiqc","bwameth_multiqc_report.html")]]
    }

}



# FIXME: add all relevant bwameth realted rules here
# if USEBWAMETH: 
#     d_targets.append('bwameth-mapping-stats')

selected_targets_default = [] 

if USEBISMARK:
    # Should we perform differential analysis?
    if config["DManalyses"]:
      selected_targets_default += ['final-report', 'diffmeth-report', 'bigwig', 'multiqc']
    else:
      selected_targets_default += ['final-report', 'bigwig','multiqc']

if USEBWAMETH:
    # Should we perform differential analysis?
    if config["DManalyses"]:
      selected_targets_default += ['final-report-bwameth', 'diffmeth-report-bwameth', 'bigwig-bwameth', 'multiqc-bwameth']
    else:
      selected_targets_default += ['final-report-bwameth', 'bigwig-bwameth','multiqc-bwameth']

# Selected output files from the above set.
selected_targets = config['execution']['target'] or selected_targets_default

# Check for availability of requested target
for target in selected_targets:
    if not target in targets.keys():
        target_desc = []
        for key in targets.keys():
            target_desc += ['{}:\t  {}'.format(key.ljust(25, ' '), targets[key]['description'])]
        
        bail("\n".join(["ERROR: Selected target '{}' is unknown.".format(target),
            "Please choose from available targets:",
            "\n"]+target_desc))

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
# Generate the multiqc report for all samples:

include: './rules/multiqc.py'


# ==========================================================================================
# Generate the final report for individual samples:

rule final_report:
    input:
        methCall_tabixfile  = os.path.join(DIR_methcall,"{tool}","tabix_{context}","{prefix}_{context}.txt.bgz"),
        methSeg_bedfile     = os.path.join(DIR_seg,"{prefix}_{context}_{tool}.meth_segments.bed"),
        template            = os.path.join(DIR_templates,"index.Rmd"),
        bigwigFile          = os.path.join(DIR_bigwig, "{prefix}.{context}_{tool}.bw")
    output:
        report        = os.path.join(DIR_final, 
                "sample_reports",
                "{prefix}_{context}_{tool}_{assembly}_final.html")
    params:
        ## absolute path to bamfiles
        sampleid                     = "{prefix}",
        source_dir                   = config['locations']['input-dir'],
        out_dir                      = OUTDIR,
        inBam                        = os.path.join(OUTDIR, DIR_sorted,"{prefix}.bam"),
        context                      = "{context}",
        assembly                     = ASSEMBLY,
        mincov                       = int(config['general']['methylation-calling']['minimum-coverage']),
        minqual                      = int(config['general']['methylation-calling']['minimum-quality']),
        TSS_plotlength               = int(config['general']['reports']['TSS_plotlength']),
        methSegPng                   = os.path.join(DIR_seg,"{prefix}_{context}_{tool}.meth_segments.png"),
        scripts_dir                  = DIR_scripts,
        refGene_bedfile              = REFGENES_BEDFILE,
        webfetch                     = WEBFETCH,
        # required for any report
        bibTexFile                   = BIBTEXPATH,
        prefix                       = "{prefix}_{context}_{tool}_{assembly}_{tool}",
        workdir                      = os.path.join(DIR_final,"sample_reports"),
        logo                         = LOGOPATH
    log:
        os.path.join(DIR_final,"sample_reports", "{prefix}_{context}_{tool}_{assembly}_final.log")
    message: fmt("Compiling final report for sample {wildcards.prefix}.")
    # run:
    #     generateReport(input, output, params, log, "")
    shell:
        nice('Rscript', ["{DIR_scripts}/generate_report.R",
                           "--reportFile={input.template}",
                           "--outFile={output.report}",
                           "--workdir={params.workdir}",
                           "--logo={params.logo}",
                           "--bibTexFile={params.bibTexFile}",
                           "--prefix={params.prefix}",
                           "--report.params='{{"+
                           ",".join([
                               '"sampleid":"{params.sampleid}"',
                               '"assembly":"{params.assembly}"',
                               '"context":"{params.context}"',
                               '"bigwigFile":"{input.bigwigFile}"',
                               '"inBam":"{params.inBam}"',
                               '"methCall_tabixfile":"{input.methCall_tabixfile}"',
                               '"methSegBed":"{input.methSeg_bedfile}"',
                               '"methSegPng":"{params.methSegPng}"',
                              
    
                               '"mincov":"{params.mincov}"',
                               '"minqual":"{params.minqual}"',
                               '"TSS_plotlength":"{params.TSS_plotlength}"',

                               '"source_dir":"{params.source_dir}"',
                               '"out_dir":"{params.out_dir}"',
                               '"scripts_dir":"{params.scripts_dir}"',
                               
                               '"refGenes_bedfile":"{params.refGene_bedfile}"',
                               '"webfetch":"{params.webfetch}"'
                           ])+"}}'",
                           "--logFile={log}"], "{log}", 
                           "echo '' ")



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
    message: fmt("Segmenting methylation profile for {wildcards.context} context of sample {wildcards.prefix}.")
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
        DIR_mapped+"{sample}_trimmed_bismark_bt2.bam"
    output:
        DIR_sorted+"{sample}_se_bt2.sorted.deduped.bam"
    params:
        threads=config['execution']['rules']['samblaster_markdup_sort']['threads'],
        memory=config['execution']['rules']['samblaster_markdup_sort']['memory'],
        tmpdir=DIR_sorted+"{sample}/"
    log:
        DIR_sorted+"{sample}_deduplication.log"
    message: fmt("Deduplicating single-end aligned reads from {input}")
    shell:
        nice("samtools", 
        ["view -h {input}"," | ", 
        tool("samblaster"),"-r",toolArgs("samblaster"),"2> {log}","|",
        tool("samtools"),"sort","-T={params.tmpdir}",
         "-o {output}", "-@ {params.threads}", 
         "-m {params.memory}", "-l 9","2> {log}",";",
         tool("samtools"),"index {output}"],("{log}"))


#-----------------------
rule deduplication_pe:
    input:
        DIR_mapped+"{sample}_1_val_1_bismark_bt2_pe.bam"
    output:
        DIR_sorted+"{sample}_1_val_1_bt2.sorted.deduped.bam"
    params:
        threads=config['execution']['rules']['samblaster_markdup_sort']['threads'],
        memory=config['execution']['rules']['samblaster_markdup_sort']['memory'],
        tmpdir=DIR_sorted+"{sample}/"
    log:
        DIR_sorted+"{sample}_deduplication.log"
    message: fmt("Deduplicating paired-end aligned reads from {input}")
    shell:
        nice("samtools", 
        ["view -h {input}"," | ", 
        tool("samblaster"),"-r",toolArgs("samblaster"),"2> {log}","|",
        tool("samtools"),"sort","-T={params.tmpdir}",
         "-o {output}", "-@ {params.threads}", 
         "-m {params.memory}", "-l 9","2> {log}",";",
         tool("samtools"),"index {output}"],("{log}"))


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
        qc     = DIR_posttrim_QC+"{sample}_trimmed_fastqc.zip"
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
        qc   = [ DIR_posttrim_QC+"{sample}_1_val_1_fastqc.zip",
                 DIR_posttrim_QC+"{sample}_2_val_2_fastqc.zip"]
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
# Merge methylation samples and perform differential analysis:

include: './rules/perform_diffmeth.py'


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
        os.path.join(OUTDIR,'bismark_genome_preparation_'+ASSEMBLY+'.log')
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
        chrom_seqlengths  = os.path.join(DIR_mapped,ASSEMBLY+"_chromlengths.csv")
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
    	DIR_posttrim_QC+"{sample}_trimmed_fastqc.zip"
    params:
        fastqc_args = config['tools']['fastqc']['args'],
        outdir = "--outdir "+DIR_posttrim_QC
    log:
   	    DIR_posttrim_QC+"{sample}_trimmed_fastqc.log"
    message: fmt("Quality checking trimmmed single-end data from sample {wildcards.sample}")
    shell:
        nice('fastqc', ["{params}", "{input}"], "{log}")

rule fastqc_after_trimming_pe:
    input:
        DIR_trimmed+"{sample}_1_val_1.fq.gz",
        DIR_trimmed+"{sample}_2_val_2.fq.gz"
    output:
    	DIR_posttrim_QC+"{sample}_1_val_1_fastqc.zip",
    	DIR_posttrim_QC+"{sample}_2_val_2_fastqc.zip"
    params:
        fastqc_args = config['tools']['fastqc']['args'],
        outdir = "--outdir "+DIR_posttrim_QC
    log:
   	    DIR_posttrim_QC+"{sample}_trimmed_fastqc.log"
    message: fmt("Quality checking trimmmed paired-end data from sample {wildcards.sample}")
    shell:
        nice('fastqc', ["{params}", "{input}"], "{log}")


# ==========================================================================================
# Trim the reads for adapter-ends and quality

rule trim_reads_se:
    input:
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
    message: fmt("Trimming raw single-end read data from sample {wildcards.sample}")
    shell:
       nice('trim-galore', ["{params}", "{input.file}"], "{log}")

rule trim_reads_pe:
    input:
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
        fmt("Trimming raw paired-end read data from sample {wildcards.sample}")
    shell:
        nice('trim-galore', ["{params}", "{input.files}"], "{log}")


# ==========================================================================================
# Perform quality control on raw data

rule fastqc_raw: #----only need one: covers BOTH pe and se cases.
    input:
        PATHIN+"{sample}.fq.gz"
    output:
        DIR_rawqc+"{sample}_fastqc.zip"
    params:
        fastqc_args = config['tools']['fastqc']['args'],
        outdir      = "--outdir "+ DIR_rawqc     # usually pass params as strings instead of wildcards.
    log:
        DIR_rawqc+"{sample}_fastqc.log"
    message: fmt("Quality checking raw read data from sample {wildcards.sample}")
    shell:
        nice('fastqc', ["{params}", "{input}"], "{log}")


