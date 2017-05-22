
def getFilenames(mylist):
  if len(mylist)==2:
    return( [splitext_fqgz(mylist[0])[0], splitext_fqgz(mylist[1])[0]] )
  if len(mylist)==1:
    return( [splitext_fqgz(mylist[0])[0]] )
  else:
    raise Exception("Sth went wrong in getFilenames())")
    
    
def getExtension(mylist):
  if len(mylist)==2:
    return( [splitext_fqgz(mylist[0])[1], splitext_fqgz(mylist[1])[1]] )
  if len(mylist)==1:
    return( [splitext_fqgz(mylist[0])[1]] )
  else:
    raise Exception("Sth went wrong in getExtension())")   
    


def list_files(PATH, files, sampleid, ext):
    if len(files) == 1:
        return [PATH+sampleid+ext] #---- single end
    elif len(files) == 2:
        return [PATH+sampleid+"_1"+ext, PATH+sampleid+"_2"+ext] #---- paired end
    else:
        raise Exception("=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")


def mylist_files_TG(PATH, files, sampleid):
    if len(files) == 1:
        return [PATH+sampleid+"_trimmed.fq.gz"] #---- single end
    elif len(files) == 2:
        return [PATH+sampleid+"_1_val_1.fq.gz", PATH+sampleid+"_2_val_2.fq.gz"] #---- paired end
    else:
        raise Exception("=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")


# def list_files_TG(PATH, files, sampleid):
#     if len(files) == 1:
#         return [PATH+files[0]+"_trimmed.fq.gz"] #---- single end
#     elif len(files) == 2:
#         return [PATH+files[0]+"_val_1.fq.gz", PATH+files[1]+"_val_2.fq.gz"] #---- paired end
#     else:
#         raise Exception("=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")
# 


def list_files_posttrim_QC(PATH, files, sampleid, ext):
    if len(files) == 1:
        return [PATH+sampleid+"_trimmed_fastqc"+ext] #---- single end
    elif len(files) == 2:
        return [PATH+sampleid+"_1_val_1_fastqc"+ext, PATH+sampleid+"_2_val_2_fastqc"+ext] #---- paired end
    else:
        raise Exception("=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")

def list_files_bismark(PATH, files, sampleid):
    return(PATH+sampleid+".bam")
    
    
# def real_list_files_bismark(PATH, files, sampleid):
#     if len(files) == 1:
#         return [
#                 #PATH+files[0]+"_trimmed_bismark_bt2_SE_report.txt", #TODO: comment it for now, just make it simpler
#                 PATH+files[0]+"_trimmed_bismark_bt2.bam"] #---- single end
#     elif len(files) == 2:
#         return [
#                 #PATH+files[0]+"_val_1_bismark_bt2_PE_report.txt",
#                 PATH+files[0]+"_val_1_bismark_bt2_pe.bam"] #---- paired end
#     else:
#         raise Exception("=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")

# 
# def list_files_dedupe(PATH, files):
#     if len(files) == 1:
#         return [PATH+files[0]+"_se.deduplicated.bam"] #---- single end
#     elif len(files) == 2:
#         return [PATH+files[0]+"_val_1.deduplicated.bam"] #---- paired end
#     else:
#         raise Exception("=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")
# 
# 
# def list_files_sortbam(PATH, files):
#     if len(files) == 1:
#         return [PATH+files[0]+"_se.deduplicated.sorted.bam"] #---- single end
#     elif len(files) == 2:
#         return [PATH+files[0]+"_val_1.deduplicated.sorted.bam"] #---- paired end
#     else:
#         raise Exception("=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")
# 
# def SEPEstr(files):
#     if len(files) == 1:
#         return  "_trimmed_bismark_bt2_SE" #---- single end
#     elif len(files) == 2:
#         return  "_val_1_bismark_bt2_PE" #---- paired end
#     else:
#         raise Exception("=== ERROR: file list is neither 1 nor 2 in length for file[0]="+files[0]+". HALTING! ===")
# 
