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


def SEPEstr(files):
    if len(files) == 1:
        return  "_trimmed_bismark_bt2_SE" #---- single end
    elif len(files) == 2:
        return  "_val_1_bismark_bt2_PE" #---- paired end
    else:
        raise Exception("=== ERROR: file list is neither 1 nor 2 in length for file[0]="+files[0]+". HALTING! ===")

def Annot(PATH, files, assembly, sampleID):
    if len(files) == 1:
      return  PATH+sampleID+"_se_bt2.deduped.sorted_"+assembly+"_annotation.nb.html" #---- single end
    elif len(files) == 2:
        return [PATH+sampleID+"_1_val_1_bt2.deduped.sorted_"+assembly+"_annotation.nb.html"] #---- paired end

def fmt(message):
    """Format the MESSAGE string."""
    return "----------  " + message + "  ----------"

def list_files_xmeth(PATH, files, sampleID):
    if len(files) == 1:
        return [ PATH+sampleID+"_se_bt2.deduped.bedGraph.gz",
                 PATH+sampleID+"_se_bt2.deduped.bismark.cov.gz",
                 PATH+sampleID+"_se_bt2.deduped.CpG_report.txt.gz"] #---- single end 
    elif len(files) == 2:
        return [PATH+sampleID+"_1_val_1_bt2.deduped.bedGraph.gz",
                PATH+sampleID+"_1_val_1_bt2.deduped.bismark.cov.gz",
                PATH+sampleID+"_1_val_1_bt2.deduped.CpG_report.txt.gz"] #---- paired end 
    else:
        raise Exception("=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")   

def bam_processing(PATH, files, sampleID):
    if len(files) == 1:
      return  PATH+sampleID+"_se_bt2.deduped.sorted_meth_calls.nb.html" #---- single end
    elif len(files) == 2:
        return [PATH+sampleID+"_1_val_1_bt2.deduped.sorted_meth_calls.nb.html"] #---- paired end

def DiffMeth(PATH, treatments):
    return [PATH+"_".join(treatments)+".sorted_diffmeth.nb.html"]

def Final(PATH, files, assembly, sampleID):
    if len(files) == 1:
      return  PATH+sampleID+"_se_bt2.deduped.sorted_"+assembly+"_final.nb.html" #---- single end
    elif len(files) == 2:
        return [PATH+sampleID+"_1_val_1_bt2.deduped.sorted_"+assembly+"_final.nb.html"] #---- paired end

        
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

def diff_meth_input(wc):
  sample = wc.prefix
  sampleid = get_fastq_name(sample)
  treatment_of_sampleid = SAMPLE_TREATMENTS_DICT[ sampleid ]

  mylist = []
  for x in DIFF_METH_TREATMENT_PAIRS:
    if treatment_of_sampleid in x:
      name_of_dir = x[0]+"_"+x[1]+".sorted_"+wc.assembly+"_annotation.diff.meth.nb.html"
      mylist.append(DIR_annot + name_of_dir)
  return(mylist)

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


def generateReport(input, output, params, log, reportSubDir):
    dumps = json.dumps(dict(params.items()),sort_keys=True,
                       separators=(",",":"), ensure_ascii=True)

    cmd =   "{RSCRIPT} {DIR_scripts}/report_functions.R"
    cmd +=  " --reportFile={input.template}"
    cmd +=  " --outFile={output.report}"
    cmd +=  " --finalReportDir=" + os.path.join(DIR_final,reportSubDir)
    cmd +=  " --report.params={dumps:q}"
    cmd +=  " --logFile={log}"
    shell(cmd, dumps)

#--- NICE gauges the computational burden, ranging from -19 to +19.
#--- The more "nice" you are, the more you allow other processes to jump ahead of you
#--- (like in traffic). Generally set to maximally nice=19 to avoid interference with other users.
def nice(cmd):
    return "nice -" + str(config['execution']['nice']) + " " + cmd
