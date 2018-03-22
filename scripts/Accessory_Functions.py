# ----------------------------------------------------------------------------- #
def java_tool(java, threads, mem, tempdir, tool_path, tool_name):
    
    # removes 1 g from java memory heap
    mem = int(float(mem[:,-1]))
    if(mem < 1):
        mem = mem/2
    else:
        mem = mem -1
    mem = mem + 'G'

    tool = ' '.join([
        java,
        '-XX:ParallelGCThreads=' + str(threads),
        '-Xmx'                   + str(mem),
        '-Djava.io.tmpdir='      + str(tempdir),
        '-jar',
        str(tool_path),
        str(tool_name)
    ])
    return tool
