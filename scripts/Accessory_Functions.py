from snakemake import shell

# ----------------------------------------------------------------------------- #
# java_head_difference is the reduction in heap allocation given to the java executible
# compared to SGE submission: if SGE submission requests 16G, java will be
# started with (16 - java_heap_difference)G
def java_tool(java, threads, mem, tempdir, tool_path, tool_name, java_heap_difference=5):

    # removes 1 g from java memory heap
    mem_size   = int(float(mem[:-1]))
    mem_suffix = mem[-1]
    if(mem_size <= java_heap_difference):
        mem_size = java_heap_difference - 1
    else:
        mem_size = mem_size - java_heap_difference

    mem_reduced = str(mem_size) + mem_suffix

    tool = ' '.join([
        java,
        '-XX:ParallelGCThreads=' + str(threads),
        '-Xmx'                   + str(mem_reduced),
        '-Xss'                   + str('16M'),
        '-Djava.io.tmpdir='      + str(tempdir),
        '-jar',
        str(tool_path),
        str(tool_name)
    ])
    return tool

# ----------------------------------------------------------------------------- # prints the command to STDERR and executes
def print_shell(command):
    print(command, file=sys.stderr)
    shell(command)
    

# ----------------------------------------------------------------------------- 
# extracts the adapter start and adapter length from the adapter hash
# sample_name is the sample name defined in sample sheet
# type : umi_barcode / cell_barcode
def adapter_params(sample_name, type):

    if not type in set(['cell_barcode','umi_barcode']):
        sys.exit('invalid barcode type')

    method = SAMPLE_SHEET.fetch_field(sample_name,'method')[0]
    adapter_params = ADAPTER_PARAMETERS[method][type]
    barcode_hash = {
        'start'  : adapter_params['base_min'],
        'length' : adapter_params['base_max'] - adapter_params['base_min'] + 1
    }
    return(barcode_hash)
# ----------------------------------------------------------------------------- # calculates the barcode length from the sample sheet
def get_adapter_size(name):
    cb_adapter   = adapter_params(name, 'cell_barcode')
    umi_adapter  = adapter_params(name, 'umi_barcode')
    adapter_size = cb_adapter['length'] + umi_adapter['length']
    return(adapter_size)
    
    