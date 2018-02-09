# ---------------------------------------------------------------------------- #
options = commandArgs(trailingOnly=TRUE)
source(file.path(options[2],'/Argument_Parser.R'))
argv = Parse_Arguments('Parse_Bowtie2log')

# ---------------------------------------------------------------------------- #
Parse_Bowtie2log = function(
    infiles,
    outfile,
    scriptdir
){

    
    library(stringr)
    source(file.path(scriptdir,'Functions_Helper.R'))
    
    # ------------------------------------------------------------------------ #
    
    message('Parsing log files ...')
    ld = list()
    for(i in 1:length(infiles)){
        
        infile = infiles[i]
        name = str_replace(basename(infile),'.bowtie2.log','')
        message(name)
        stat = MappingStats_Bowtie2(infile)
        stat$sample_name = name
        ld[[name]] = stat
        
    }
    dd = do.call(rbind, ld)
    saveRDS(dd, file = outfile)
}

# ---------------------------------------------------------------------------- #
# function call
Parse_Bowtie2log(
    infiles     = argv$params[['logfiles']],
    outfile     = argv$output[['outfile']],
    scriptdir   = argv$params[['scriptdir']]
)
