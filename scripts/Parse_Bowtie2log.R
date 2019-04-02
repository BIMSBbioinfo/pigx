# ---------------------------------------------------------------------------- #
options = commandArgs(trailingOnly=TRUE)
source(file.path(options[2],'/Argument_Parser.R'))
argv = Parse_Arguments('Parse_Bowtie2log')

# ---------------------------------------------------------------------------- #
Parse_Bowtie2log = function(
    path_log,
    outfile,
    genome_types,
    scriptdir
){
    library(stringr)
    source(file.path(scriptdir,'Functions_Helper.R'))

    # ------------------------------------------------------------------------ #

    message('Parsing log files ...')
    infiles = list.files(path_log, full.names=TRUE, recursive = TRUE ,pattern='bowtie2')
    infiles = infiles[!str_detect(infiles,'build')]

    ld = list()
    for(genome_type in genome_types){

        infiles_genome = infiles[str_detect(infiles, genome_type)]

        # checks whether any file was mapped to the corresponding genome
        if(length(infiles_genome) > 0){
          for(i in 1:length(infiles_genome)){

              infile = infiles_genome[i]
              name = str_replace(basename(infile),'.bowtie2.+.log','')
              message(name)
              stat = MappingStats_Bowtie2(infile)
              stat$sample_name = name
              stat$genome_type = genome_type
              ld = c(ld, list(stat))

          }
        }
    }
    dd = do.call(rbind, ld)
    saveRDS(dd, file = outfile)

    # checks whether the parsing succeeded
    if(file.size(outfile) < 100){
      stop('Bowtie2 log parse is empty')
    }else{
      return(1)
    }
}

# ---------------------------------------------------------------------------- #
# function call
Parse_Bowtie2log(
    path_log     = argv$params[['path_log']],
    outfile      = argv$output[['outfile']],
    scriptdir    = argv$params[['scriptdir']],
    genome_types = argv$params[['genome_types']]
)
