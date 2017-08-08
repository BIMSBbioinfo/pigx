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

def Annot(PATH, files, assembly):
    if len(files) == 1:
      return  PATH+sampleID+"_se_bt2.deduped.sorted_"+assembly+"_annotation.nb.html" #---- single end
    elif len(files) == 2:
        return [PATH+sampleID+"_1_val_1_bt2.deduped.sorted_"+assembly+"_annotation.nb.html"] #---- paired end

def fmt(message):
    """Format the MESSAGE string."""
    return "----------  " + message + "  ----------"
