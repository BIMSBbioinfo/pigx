def list_files(PATH, files, ext):
    if len(files) == 1:
        return [PATH+files[0]+ext] #---- single end
    elif len(files) == 2:
        return [PATH+files[0]+ext, PATH+files[1]+ext] #---- paired end
    else:
        raise Exception("=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")


def list_files_TG(PATH, files):
    if len(files) == 1:
        return [PATH+files[0]+"_trimmed.fq.gz"] #---- single end
    elif len(files) == 2:
        return [PATH+files[0]+"_val_1.fq.gz", PATH+files[1]+"_val_2.fq.gz"] #---- paired end
    else:
        raise Exception("=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")


def list_files_posttrim_QC(PATH, files, ext):
    if len(files) == 1:
        return [PATH+files[0]+"_trimmed_fastqc"+ext] #---- single end
    elif len(files) == 2:
        return [PATH+files[0]+"_val_1_fastqc"+ext, PATH+files[1]+"_val_2_fastqc"+ext] #---- paired end
    else:
        raise Exception("=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")

def list_files_bismark(PATH, files):
    if len(files) == 1:
        return [PATH+files[0]+"_trimmed_bismark_bt2_SE_report.txt",
                PATH+files[0]+"_trimmed_bismark_bt2.bam", 
                PATH+files[0]+"_trimmed_bismark_bt2.nucleotide_stats.txt" ] #---- single end
    elif len(files) == 2:
        return [PATH+files[0]+"_val_1_bismark_bt2_PE_report.txt",
                PATH+files[0]+"_val_1_bismark_bt2_pe.bam",
                PATH+files[0]+"_val_1_bismark_bt2_pe.nucleotide_stats.txt"] #---- paired end
    else:
        raise Exception("=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")


def list_files_dedupe(PATH, files):
    if len(files) == 1:
        return [PATH+files[0]+"_se.deduplicated.bam"] #---- single end
    elif len(files) == 2:
        return [PATH+files[0]+"_val_1.deduplicated.bam"] #---- paired end
    else:
        raise Exception("=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")


def list_files_sortbam(PATH, files):
    if len(files) == 1:
        return [PATH+files[0]+"_se.deduplicated.sorted.bam"] #---- single end
    elif len(files) == 2:
        return [PATH+files[0]+"_val_1.deduplicated.sorted.bam"] #---- paired end
    else:
        raise Exception("=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")


def list_files_deconv(PATH, files):
    if len(files) == 1:
        return [ PATH+files[0]+"_se.deconv_out.RData"] #---- single end
    elif len(files) == 2:
        return [PATH+files[0]+"_val_1.deconv_out.RData"] #---- paired end
    else:
        raise Exception("=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")

def SEPEstr(files):
    if len(files) == 1:
        return  "_trimmed_bismark_bt2_SE" #---- single end
    elif len(files) == 2:
        return  "_val_1_bismark_bt2_PE" #---- paired end
    else:
        raise Exception("=== ERROR: file list is neither 1 nor 2 in length for file[0]="+files[0]+". HALTING! ===")

