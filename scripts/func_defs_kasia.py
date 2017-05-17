def list_files(PATH, files, ext):
    if len(files) == 1:
        return [PATH+files[0]+ext] #---- single end
    elif len(files) == 2:
        return [PATH+files[0]+ext, PATH+files[1]+ext] #---- paired end
    else:
        raise Exception("=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")

def list_files_TG(PATH, files, ext):
    if len(files) == 1:
        return [PATH+files[0]+"_trimmed"+ext] #---- single end
    elif len(files) == 2:
        return [PATH+files[0]+"_val_1"+ext, PATH+files[1]+"_val_2"+ext] #---- paired end
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
                PATH+files[0]+"_trimmed_bismark_bt2.bam" ] #---- single end
    elif len(files) == 2:
        return [PATH+files[0]+"_val_1_bismark_bt2_PE_report.txt",
                PATH+files[0]+"_val_1_bismark_bt2_pe.bam"] #---- paired end
    else:
        raise Exception("=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")


def list_files_dedupe(PATH, files):
    if len(files) == 1:
        return [PATH+files[0]+"_trimmed_bismark_bt2.deduplicated.bam",
                PATH+files[0]+"_trimmed_bismark_bt2.deduplication_report.txt"] #---- single end
    elif len(files) == 2:
        return [PATH+files[0]+"_val_1_bismark_bt2_pe.deduplicated.bam",
                PATH+files[0]+"_val_1_bismark_bt2_pe.deduplication_report.txt"] #---- paired end
    else:
        raise Exception("=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")

def list_files_sortbam(PATH, files):
    if len(files) == 1:
        return [PATH+files[0]+"_trimmed_bismark_bt2.deduplicated.sorted.bam"] #---- single end
    elif len(files) == 2:
        return [PATH+files[0]+"_val_1_bismark_bt2_pe.deduplicated.sorted.bam"] #---- paired end
    else:
        raise Exception("=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")


def list_files_indexbam(PATH, files):
    if len(files) == 1:
        return [PATH+files[0]+"_trimmed_bismark_bt2.deduplicated.sorted.bam.bai"] #---- single end
    elif len(files) == 2:
        return [PATH+files[0]+"_val_1_bismark_bt2_pe.deduplicated.sorted.bam.bai"] #---- paired end
    else:
        raise Exception("=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")


