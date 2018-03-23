# ----------------------------------------------------------------------------- #
# java_head_difference is the reduction in heap allocation given to the java executible
# compared to SGE submission: if SGE submission requests 16G, java will be 
# started with (16 - java_heap_difference)G
def java_tool(java, threads, mem, tempdir, tool_path, tool_name, java_heap_difference=1):
    
    # removes 1 g from java memory heap
    mem_size   = int(float(mem[:-1]))
    mem_suffix = mem[-1]
    if(mem_size < java_heap_difference):
        mem_size = mem_size/2
    else:
        mem_size = mem_size - java_heap_difference
    
    mem_reduced = str(mem) + mem_suffix

    tool = ' '.join([
        java,
        '-XX:ParallelGCThreads=' + str(threads),
        '-Xmx'                   + str(mem_reduced),
        '-Djava.io.tmpdir='      + str(tempdir),
        '-jar',
        str(tool_path),
        str(tool_name)
    ])
    return tool
