# ----------------------------------------------------------------------------- #
def java_tool(java, threads, mem, tempdir, droptools, name):
    tool = ' '.join([
        java,
        '-XX:ParallelGCThreads=' + threads,
        '-Xmx$'                  + mem,
        '-Djava.io.tmpdir='      + tempdir,
        droptools,
        name
    ])
    return tool
