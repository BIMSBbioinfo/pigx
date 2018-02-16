# ----------------------------------------------------------------------------- #
def java_tool(java, threads, mem, tempdir, tool_path, tool_name):
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
