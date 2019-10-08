# PiGx BSseq Pipeline.
#
# Copyright © 2017 Bren Osberg <brendan.osberg@mdc-berlin.de>
# Copyright © 2017, 2018 Alexander Gosdschan <alexander.gosdschan@mdc-berlin.de>
# Copyright © 2017 Katarzyna Wreczycka <katwre@gmail.com>
# Copyright © 2017, 2018 Ricardo Wurmus <ricardo.wurmus@mdc-berlin.de>
#
# This file is part of the PiGx BSseq Pipeline.
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
from glob import glob

def dedupe_tag(protocol):
    if protocol.upper() == "WGBS":
        return ".deduped"
    elif protocol.upper() == "RRBS":
        return ""
    else:
        raise Exception("=== ERROR: unexpected protocol ===")

def files_for_sample(proc):
    return [expand(proc(config['SAMPLES'][sample]['files'],
                        config['SAMPLES'][sample]['SampleID'],
                        config['SAMPLES'][sample]['Protocol']))
            for sample in config['SAMPLES']]


def list_files_rawQC(files, sampleID, protocol):
    PATH = DIR_rawqc
    if len(files) == 1:
        return [PATH+sampleID+"_fastqc.html"]  # ---- single end
    elif len(files) == 2:
        # ---- paired end
        return [PATH+sampleID+"_1_fastqc.html", PATH+sampleID+"_2_fastqc.html"]
    else:
        raise Exception(
            "=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")


def list_files_TG(files, sampleID, protocol):
    PATH = DIR_trimmed
    if len(files) == 1:
        return [PATH+sampleID+"_trimmed.fq.gz"]  # ---- single end
    elif len(files) == 2:
        # ---- paired end
        return [PATH+sampleID+"_1_val_1.fq.gz", PATH+sampleID+"_2_val_2.fq.gz"]
    else:
        raise Exception(
            "=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")


def list_files_posttrim_QC(files, sampleID, protocol):
    PATH = DIR_posttrim_QC
    if len(files) == 1:
        return [PATH+sampleID+"_trimmed_fastqc.html"]  # ---- single end
    elif len(files) == 2:
        # ---- paired end
        return [PATH+sampleID+"_1_val_1_fastqc.html", PATH+sampleID+"_2_val_2_fastqc.html"]
    else:
        raise Exception(
            "=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")


def list_files_bismark(files, sampleID, protocol):
    PATH = DIR_mapped
    if len(files) == 1:
        return [PATH+sampleID+"_trimmed_bismark_bt2_SE_report.txt",
                PATH+sampleID+"_trimmed_bismark_bt2.bam"]  # ---- single end
    elif len(files) == 2:
        return [PATH+sampleID+"_1_val_1_bismark_bt2_PE_report.txt",
                PATH+sampleID+"_1_val_1_bismark_bt2_pe.bam"]  # ---- paired end
    else:
        raise Exception(
            "=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")


def list_files_bwameth(files, sampleID, protocol):
    PATH = DIR_mapped
    if len(files) == 1:
        return [PATH+sampleID+".bwameth.bam"]  # ---- single end
    elif len(files) == 2:
        return [PATH+sampleID+".bwameth.bam"]  # ---- paired end
    else:
        raise Exception(
            "=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")


def list_files_dedupe(files, sampleID, protocol):
    PATH = DIR_sorted
    if len(files) == 1:
        # ---- single end
        return [PATH+sampleID+"_se_bt2.sorted" + dedupe_tag(protocol) + ".bam"]
    elif len(files) == 2:
        # ---- paired end
        return [PATH+sampleID+"_1_val_1_bt2.sorted" + dedupe_tag(protocol) + ".bam"]
    else:
        raise Exception(
            "=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")


def list_files_markdup(files, sampleID, protocol):
    PATH = DIR_mapped
    if len(files) == 1:
        return [PATH+sampleID+".bwameth.sorted.markdup.bam"]  # ---- single end
    elif len(files) == 2:
        return [PATH+sampleID+".bwameth.sorted.markdup.bam"]  # ---- paired end
    else:
        raise Exception(
            "=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")


def list_files_bwamethMappingStats(files, sampleID, protocol):
    PATH = DIR_mapped
    if len(files) == 1:
        return [PATH+sampleID+".bwameth.sorted.markdup.idxstats.txt",
                PATH+sampleID+".bwameth.sorted.markdup.stats.txt",
                PATH+sampleID+".bwameth.sorted.markdup.flagstat.txt"]  # ---- single end
    elif len(files) == 2:
        return [PATH+sampleID+".bwameth.sorted.markdup.idxstats.txt",
                PATH+sampleID+".bwameth.sorted.markdup.stats.txt",
                PATH+sampleID+".bwameth.sorted.markdup.flagstat.txt"]  # ---- paired end
    else:
        raise Exception(
            "=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")


def list_files_sortbam(files, sampleID, protocol):
    PATH = DIR_sorted
    if len(files) == 1:
        return [PATH+sampleID+"_se_bt2.sorted.bam"]  # ---- single end
    elif len(files) == 2:
        return [PATH+sampleID+"_1_val_1_bt2.sorted.bam"]  # ---- paired end
    else:
        raise Exception(
            "=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")


def bam_processing(files, sampleID, protocol):
    PATH = DIR_methcall
    if len(files) == 1:
        return  PATH+sampleID+"_se_bt2.sorted" + dedupe_tag(protocol) + "_methylRaw.RDS" #---- single end
    elif len(files) == 2:
        return [PATH+sampleID+"_1_val_1_bt2.sorted" + dedupe_tag(protocol) + "_methylRaw.RDS"] #---- paired end

def list_files_methyldackel_extract(files, sampleID, protocol):
    PATH = DIR_methcall + "methylDackel/"
    if len(files) == 1:
        # ---- single end
        return [PATH+sampleID+dedupe_tag(protocol)+"_methyldackel_CpG.methylKit",
                PATH+sampleID+dedupe_tag(protocol)+"_methyldackel_CHG.methylKit",
                PATH+sampleID+dedupe_tag(protocol)+"_methyldackel_CHH.methylKit",
                PATH+sampleID+"_mbias_methyldackel.txt",
                PATH+sampleID+"_mbias_OB.svg",
                PATH+sampleID+"_mbias_OT.svg",
                PATH+sampleID+"_methyldackel.cytosine_report.txt"
]
    elif len(files) == 2:
        return [PATH+sampleID+dedupe_tag(protocol)+"_methyldackel_CpG.methylKit",
                PATH+sampleID+dedupe_tag(protocol)+"_methyldackel_CHG.methylKit",
                PATH+sampleID+dedupe_tag(protocol)+"_methyldackel_CHH.methylKit",
                PATH+sampleID+"_mbias_methyldackel.txt",
                PATH+sampleID+"_mbias_OB.svg",
                PATH+sampleID+"_mbias_OT.svg",
                PATH+sampleID+"_methyldackel.cytosine_report.txt"
                ]  # ---- paired end
    else:
        raise Exception(
            "=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")

def list_files_maketabix_methyldackel(files, sampleID, protocol):
    PATH = DIR_methcall + "methylDackel/"
    if len(files) == 1:
        # ---- single end
        return [PATH+"tabix_CpG/"+sampleID+dedupe_tag(protocol)+"_CpG.txt.bgz",
                PATH+"tabix_CpG/"+sampleID+dedupe_tag(protocol)+"_CpG.txt.bgz.tbi"
                ]
    elif len(files) == 2:
        # ---- paired end
        return [PATH+"tabix_CpG/"+sampleID+dedupe_tag(protocol)+"_CpG.txt.bgz",
                PATH+"tabix_CpG/"+sampleID+dedupe_tag(protocol)+"_CpG.txt.bgz.tbi"
                ]      else:
        raise Exception(
            "=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")

def bigwig_exporting(files, sampleID, protocol):
    PATH = os.path.join(config['locations']['output-dir'], DIR_bigwig )
    if len(files) == 1:
        return  PATH+sampleID+"_se_bt2.sorted" + dedupe_tag(protocol) + ".bw" #---- single end
    elif len(files) == 2:
        return [PATH+sampleID+"_1_val_1_bt2.sorted" + dedupe_tag(protocol) + ".bw"] #---- paired end

def methSeg(files, sampleID, protocol):
    PATH = DIR_seg
    if len(files) == 1:
        return  PATH+sampleID+"_se_bt2.sorted" + dedupe_tag(protocol) + "_meth_segments_gr.RDS" #---- single end
    elif len(files) == 2:
        return [PATH+sampleID+"_1_val_1_bt2.sorted" + dedupe_tag(protocol) + "_meth_segments_gr.RDS"] #---- paired end
    else:
        raise Exception("=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")

def list_final_reports(files, sampleID, protocol):
    PATH = DIR_final
    if len(files) == 1:
        return  PATH+sampleID+"_se_bt2.sorted" + dedupe_tag(protocol) + "_"+ASSEMBLY+"_final.html" #---- single end
    elif len(files) == 2:
        return [PATH+sampleID+"_1_val_1_bt2.sorted" + dedupe_tag(protocol) + "_"+ASSEMBLY+"_final.html"] #---- paired end
    else:
        raise Exception("=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")

def fmt(message):
    """Format the MESSAGE string."""
    return "----------  " + message + "  ----------"

def get_fastq_name(full_name):
    # single end
    find_se_inx=full_name.find('_se_bt2')
    # paired-end
    find_pe_inx=full_name.find('_1_val_1_bt2')

    if(find_se_inx>=0):
      output=full_name[:find_se_inx]
    elif(find_pe_inx>=0):
     output=full_name[:find_pe_inx]
    else:
     bail("Unable to infer sample fastq name; cannot find trimming string in filename. \nHave the files been trimmed for adapter sequence and quality?")

    return(output)

def get_sampleids_from_treatment(treatment):
  treatment = treatment.replace(".deduped", "")

  SAMPLE_IDS = list(config["SAMPLES"].keys())
  SAMPLE_TREATMENTS = [config["SAMPLES"][s]["Treatment"] for s in SAMPLE_IDS]

  treatments = treatment.split("_")
  sampleids_list = []
  for t in treatments:
    sampleids = [SAMPLE_IDS[i] for i, x in enumerate(SAMPLE_TREATMENTS) if x == t]
    sampleids_list.append(sampleids)

  sampleids_list = list(sum(sampleids_list, []))
  return(sampleids_list)

def makeDiffMethPath(DIR_diffmeth, suffix, wc):
    return DIR_diffmeth + str(wc.treatment).replace('vs', '_') + dedupe_tag(config["SAMPLES"][get_sampleids_from_treatment(wc.treatment[0])[0]]['Protocol']) + '_' + suffix

# For only CpG context
def diffmeth_input_function(wc):
  treatments = wc.treatment
  treatments = treatments.replace(".deduped", "")
  sampleids  = get_sampleids_from_treatment(treatments)

  inputfiles = []
  for sampleid in sampleids:
    fqname = config["SAMPLES"][sampleid]['fastq_name']
    protocol = config["SAMPLES"][sampleid]['Protocol']
    if len(fqname)==1:
      inputfile=[os.path.join(DIR_methcall,sampleid+"_se_bt2.sorted" + dedupe_tag(protocol) + "_methylRaw.RDS")]
    elif len(fqname)==2:
      inputfile=[os.path.join(DIR_methcall,sampleid+"_1_val_1_bt2.sorted" + dedupe_tag(protocol) + "_methylRaw.RDS")]
    inputfiles.append(inputfile)

  inputfiles = list(sum(inputfiles, []))
  return(inputfiles)

def tool(name):
    return config['tools'][name]['executable']

def toolArgs(name):
    if 'args' in config['tools'][name]:
        return config['tools'][name]['args']
    else:
        return ""

# Generate a command line string that can be passed to snakemake's
# "shell".  The string is prefixed with an invocation of "nice".
def nice(cmd, args, log=None):
    executable = tool(cmd)
    line = ["nice -" + str(config['execution']['nice']), executable] + [toolArgs(cmd)] + args
    if log:
        line.append("> {} 2>&1".format(log))
    return " ".join(line)

# custome function for submitting shell command to generate final reports:
def generateReport(input, output, params, log, reportSubDir):
    dumps = json.dumps(dict(params.items()),sort_keys=True,
                       separators=(",",":"), ensure_ascii=True)

    cmd = nice('Rscript', ["{DIR_scripts}/generate_report.R",
                           "--scriptsDir=" + DIR_scripts,
                           "--reportFile={input.template}",
                           "--outFile={output.report}",
                           "--finalReportDir=" + os.path.join(DIR_final,reportSubDir),
                           "--report.params={dumps:q}",
                           "--logFile={log}"])
    print("==== The present shell command being submitted by the generateReport function is as follows: ")
    print(cmd)
    print("==== ... and the dumps string is : ")
    print(dumps)
    print("==== End of dump string.")
    print("==== The arguments list for the rendering is provided in the *html.RenderArgs.rds file ")

    shell(cmd, dumps)

# abandone current execution with a helpful error message:
def bail(msg):
    """Print the error message to stderr and exit."""
    print(msg, file=sys.stderr)
    exit(1)

# check for common input/configuration errors:
def validate_config(config):
    # Check that all locations exist
    for loc in config['locations']:
        if (not loc == 'output-dir') and (not os.path.isdir(config['locations'][loc]) ) :
            bail("ERROR: The following necessary directory does not exist: {} ({})".format(config['locations'][loc], loc))

    # Check that all of the requested differential methylation
    # treatment values are found in the sample sheet.
    treatments = set([ config["SAMPLES"][sample]["Treatment"] for sample in config["SAMPLES"] ])
    for pair in config['general']['differential-methylation']['treatment-groups']:
        for group in pair:
            if ( (group not in treatments) or (not (str.isdigit(group))) )  :
                bail("ERROR: Invalid treatment group '{}' in pair '{}'".format(group, pair))

    # Check for a genome fasta file
    fasta = glob(os.path.join(config['locations']['genome-dir'], '*.fasta'))
    fa    = glob(os.path.join(config['locations']['genome-dir'], '*.fa'))

    # Check if we have permission to write to the reference-genome directory ourselves
    # if not, then check if the ref genome has already been converted
    if (not os.access(config['locations']['genome-dir'], os.W_OK) and
        not os.path.isdir(os.path.join(config['locations']['genome-dir'], 'Bisulfite_Genome'))):
        bail("ERROR: reference genome has not been bisulfite-converted, and PiGx does not have permission to write to that directory. Please either (a) provide Bisulfite_Genome conversion directory yourself, or (b) enable write permission in {} so that PiGx can do so on its own.".format(config['locations']['genome-dir']))

    if not len(fasta) + len(fa) == 1 :
        bail("ERROR: Missing (or ambiguous) reference genome: The number of files ending in either '.fasta' or '.fa' in the following genome directory does not equal one: {}".format(config['locations']['genome-dir']))

