# PiGx BSseq Pipeline.
#
# Copyright © 2017 Bren Osberg <b.osberg@tum.de>
# Copyright © 2017 Alexander Gosdschan <alexander.gosdschan@mdc-berlin.de>
# Copyright © 2017 Katarzyna Wreczycka <katwre@gmail.com>
# Copyright © 2017 Ricardo Wurmus <ricardo.wurmus@mdc-berlin.de>
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

def files_for_sample(directory, proc, samples):
    return [ expand(proc(directory, samples[sample]['files'], samples[sample]['SampleID'])) for sample in samples ]

def list_files_rawQC(PATH, files, sampleID):
    if len(files) == 1:
        return [PATH+sampleID+"_fastqc.html"] #---- single end
    elif len(files) == 2:
        return [PATH+sampleID+"_1_fastqc.html", PATH+sampleID+"_2_fastqc.html"] #---- paired end
    else:
        raise Exception("=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")


def list_files_TG(PATH, files, sampleID):
    if len(files) == 1:
        return [PATH+sampleID+"_trimmed.fq.gz"] #---- single end
    elif len(files) == 2:
        return [PATH+sampleID+"_1_val_1.fq.gz", PATH+sampleID+"_2_val_2.fq.gz"] #---- paired end
    else:
        raise Exception("=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")


def list_files_posttrim_QC(PATH, files, sampleID):
    if len(files) == 1:
        return [PATH+sampleID+"_trimmed_fastqc.html" ] #---- single end
    elif len(files) == 2:
        return [PATH+sampleID+"_1_val_1_fastqc.html", PATH+sampleID+"_2_val_2_fastqc.html"] #---- paired end
    else:
        raise Exception("=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")

def list_files_bismark(PATH, files, sampleID):
    if len(files) == 1:
        return [PATH+sampleID+"_trimmed_bismark_bt2_SE_report.txt",
                PATH+sampleID+"_trimmed_bismark_bt2.bam"] #---- single end
    elif len(files) == 2:
        return [PATH+sampleID+"_1_val_1_bismark_bt2_PE_report.txt",
                PATH+sampleID+"_1_val_1_bismark_bt2_pe.bam"] #---- paired end
    else:
        raise Exception("=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")

def list_files_dedupe(PATH, files, sampleID):
    if len(files) == 1:
        return [PATH+sampleID+"_se_bt2.deduped.bam"] #---- single end
    elif len(files) == 2:
        return [PATH+sampleID+"_1_val_1_bt2.deduped.bam"] #---- paired end
    else:
        raise Exception("=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")

def list_files_sortbam(PATH, files, sampleID):
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

        
