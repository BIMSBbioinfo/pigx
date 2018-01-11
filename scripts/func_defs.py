# PiGx BSseq Pipeline.
#
# Copyright © 2017 Bren Osberg <b.osberg@tum.de>
# Copyright © 2017 Alexander Gosdschan <alexander.gosdschan@mdc-berlin.de>
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

def files_for_sample(proc):
    return [ expand(proc(config['SAMPLES'][sample]['files'], config['SAMPLES'][sample]['SampleID'])) for sample in config['SAMPLES'] ]

def list_files_rawQC(files, sampleID):
    PATH = DIR_rawqc
    if len(files) == 1:
        return [PATH+sampleID+"_fastqc.html"] #---- single end
    elif len(files) == 2:
        return [PATH+sampleID+"_1_fastqc.html", PATH+sampleID+"_2_fastqc.html"] #---- paired end
    else:
        raise Exception("=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")


def list_files_TG(files, sampleID):
    PATH = DIR_trimmed
    if len(files) == 1:
        return [PATH+sampleID+"_trimmed.fq.gz"] #---- single end
    elif len(files) == 2:
        return [PATH+sampleID+"_1_val_1.fq.gz", PATH+sampleID+"_2_val_2.fq.gz"] #---- paired end
    else:
        raise Exception("=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")


def list_files_posttrim_QC(files, sampleID):
    PATH = DIR_posttrim_QC
    if len(files) == 1:
        return [PATH+sampleID+"_trimmed_fastqc.html" ] #---- single end
    elif len(files) == 2:
        return [PATH+sampleID+"_1_val_1_fastqc.html", PATH+sampleID+"_2_val_2_fastqc.html"] #---- paired end
    else:
        raise Exception("=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")

def list_files_bismark(files, sampleID):
    PATH = DIR_mapped
    if len(files) == 1:
        return [PATH+sampleID+"_trimmed_bismark_bt2_SE_report.txt",
                PATH+sampleID+"_trimmed_bismark_bt2.bam"] #---- single end
    elif len(files) == 2:
        return [PATH+sampleID+"_1_val_1_bismark_bt2_PE_report.txt",
                PATH+sampleID+"_1_val_1_bismark_bt2_pe.bam"] #---- paired end
    else:
        raise Exception("=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")

def list_files_dedupe(files, sampleID):
    PATH = DIR_deduped
    if len(files) == 1:
        return [PATH+sampleID+"_se_bt2.deduped.bam"] #---- single end
    elif len(files) == 2:
        return [PATH+sampleID+"_1_val_1_bt2.deduped.bam"] #---- paired end
    else:
        raise Exception("=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")

def list_files_sortbam(files, sampleID):
    PATH = DIR_sorted
    if len(files) == 1:
        return [PATH+sampleID+"_se_bt2.deduped.sorted.bam"] #---- single end
    elif len(files) == 2:
        return [PATH+sampleID+"_1_val_1_bt2.deduped.sorted.bam"] #---- paired end
    else:
        raise Exception("=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")

def bam_processing(files, sampleID):
    PATH = DIR_methcall
    if len(files) == 1:
        return  PATH+sampleID+"_se_bt2.deduped.sorted_methylRaw.RDS" #---- single end
    elif len(files) == 2:
        return [PATH+sampleID+"_1_val_1_bt2.deduped.sorted_methylRaw.RDS"] #---- paired end

def bigwig_exporting(files, sampleID):
    PATH = DIR_bigwig
    if len(files) == 1:
        return  PATH+sampleID+"_se.bw" #---- single end
    elif len(files) == 2:
        return [PATH+sampleID+"_pe.bw"] #---- paired end
        
def methSeg(files, sampleID):
    PATH = DIR_seg
    if len(files) == 1:
        return  PATH+sampleID+"_se_bt2.deduped.sorted_meth_segments_gr.RDS" #---- single end
    elif len(files) == 2:
        return [PATH+sampleID+"_1_val_1_bt2.deduped.sorted_meth_segments_gr.RDS"] #---- paired end
        
def methSegAnnot(files, sampleID):
    PATH = DIR_seg
    if len(files) == 1:
        return  PATH+sampleID+"_se_bt2.deduped.sorted_"+ASSEMBLY+"_annotation.nb.html" #---- single end
    elif len(files) == 2:
        return [PATH+sampleID+"_1_val_1_bt2.deduped.sorted_"+ASSEMBLY+"_annotation.nb.html"] #---- paired end

def list_final_reports(files, sampleID):
    PATH = DIR_final
    if len(files) == 1:
        return  PATH+sampleID+"_se_bt2.deduped.sorted_"+assembly+"_final.html" #---- single end
    elif len(files) == 2:
        return [PATH+sampleID+"_1_val_1_bt2.deduped.sorted_"+assembly+"_final.html"] #---- paired end



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
     print("Sth went wrong")

    return(output)

SAMPLE_IDS = list(config["SAMPLES"].keys())
SAMPLE_TREATMENTS = [config["SAMPLES"][s]["Treatment"] for s in SAMPLE_IDS]

def diff_meth_input(wc):
  sample = wc.prefix
  sampleid = get_fastq_name(sample)
  sample_treatments_dict = dict(zip(SAMPLE_IDS, SAMPLE_TREATMENTS))
  treatment_of_sampleid = sample_treatments_dict[ sampleid ]

  mylist = []
  for x in config['DIFF_METH']:
    if treatment_of_sampleid in x:
      name_of_dir = x[0]+"_"+x[1]+".sorted_"+wc.assembly+"_annotation.diff.meth.nb.html"
      mylist.append(DIR_annot + name_of_dir)
  return(mylist)
  
def finalReportDiffMeth_input(prefix):
  sampleid = get_fastq_name(prefix)
  treatment_of_sampleid = SAMPLE_TREATMENTS_DICT[ sampleid ]
  treatments = ["_".join(pair) for pair in DIFF_METH_TREATMENT_PAIRS if treatment_of_sampleid in pair]
  outList = []
  if treatments: 
      outList  = [ "{}{}.sorted_{}.RDS".format(DIR_diffmeth,treat,type) for type in ["diffmeth","diffmethhyper","diffmethhypo"] for treat in treatments]
      outList += [ "{}{}.sorted_diffmeth.bed".format(DIR_diffmeth,treat,type) for treat in treatments ]
  
  return  outList

def get_sampleids_from_treatment(treatment):
  treatments = treatment.split("_")
  sampleids_list = []
  for t in treatments:
    sampleids = [SAMPLE_IDS[i] for i, x in enumerate(SAMPLE_TREATMENTS) if x == t]
    sampleids_list.append(sampleids)

  sampleids_list = list(sum(sampleids_list, []))
  return(sampleids_list)

# For only CpG context
def diffmeth_input_function(wc):
  treatments = wc.treatment
  sampleids = get_sampleids_from_treatment(treatments)

  inputfiles = []
  for sampleid in sampleids:
    fqname = config["SAMPLES"][sampleid]['fastq_name']
    if len(fqname)==1:
      inputfile=[os.path.join(DIR_methcall,sampleid+"_se_bt2.deduped.sorted_CpG.txt")]
    elif len(fqname)==2:
      inputfile=[os.path.join(DIR_methcall,sampleid+"_1_val_1_bt2.deduped.sorted_CpG.txt")]
    inputfiles.append(inputfile)

  inputfiles = list(sum(inputfiles, []))
  return(inputfiles)

def tool(name):
    return config['tools'][name]['executable']

# Generate a command line string that can be passed to snakemake's
# "shell".  The string is prefixed with an invocation of "nice".
def nice(cmd, args, log=None):
    executable = tool(cmd)
    line = ["nice -" + str(config['execution']['nice']), executable] + args
    if log:
        line.append("> {} 2>&1".format(log))
    return " ".join(line)

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
    shell(cmd, dumps)
