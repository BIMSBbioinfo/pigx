# ---------------------------------------------------------------------------- #
Make_UCSC_HUB = function(path_hub, genome_name, hub, paths){

    library(stringr)
    
    join_string = function(l, prefix=NULL){
        if(! is.null(prefix)){
            string = paste(prefix, names(l), l)
        }else{
            string = paste(names(l), l)
        }
        paste0(paste(string, collapse='\n'),'\n\n')
    }
    
    # checks for hub list
    # ------------------------------------------------------------------------ #
    # genomes.txt
    message('genomes.txt')
    file_genome = file.path(path_hub,'genomes.txt')
    genome = list(
        'genome'   = genome_name,
        'trackDb' = file.path(genome_name,'trackDb.txt'))
    cat(join_string(genome), file = file_genome)
    
    # hub.txt
    message('hub.txt')
    file_hub = file.path(path_hub,'hub.txt')
    hub_desc = list('hub'         = hub$name,
                    'shortLabel'  = hub$shortLabel,
                    'longLabel'   = hub$longLabel,
                    'genomesFile' = 'genomes.txt',
                    'email'       = hub$email,
                    'description' = hub$descriptionUrl)
    cat(join_string(hub_desc), file = file_hub)
    
    # ------------------------------------------------------------------------ #
    # track db
    message('trackDb.txt')
    path_files = file.path(path_hub, genome_name)
        dir.create(path_files)
    
    file_trackDb = file.path(path_files,'trackDb.txt')
    if(file.exists(file_trackDb))
        file.remove(file_trackDb)
    
    # ------------------------------------------------------------------------ #
    message('Supertracks ...')
    for(i in 1:length(hub$super_tracks)){
    
        super_name = names(hub$super_tracks)[i]
        message(super_name)
        long_label = paste(super_name, 'Long')
        
        tracks = hub$super_tracks[[super_name]]
        if(!is.null(tracks$long_label))
            long_label = tracks$long_label  
        
        header = list(
        'track'      = super_name, 
        'superTrack' = 'on show',
        'shortLabel' = super_name,
        'longLabel'  = long_label,
        'visibility' = 'full',
         priority    = i)
        cat(join_string(header), file = file_trackDb, append=TRUE)
        
        track_names = setdiff(names(tracks),'long_label')
        for(tn in 1:length(track_names)){
            
            track = track_names[tn]
            track_name = tracks[[track]]$name
            track_type = tracks[[track]]$type
            inpath = paths[[track_type]]$path
            suffix = paths[[track_type]]$suffix
            type   = paths[[track_type]]$type
            
            infile  = file.path(inpath,track_name, paste(track_name, suffix, sep='.'))
            if(! file.exists(infile))
                stop(paste(infile, 'does not exist'))
            outfile = file.path(path_files, basename(infile))
            file.copy(infile, outfile)
            if(! file.exists(outfile))
                stop(paste(outfile, 'does not exist'))
            
            track_description = list(
            'track'      = track_name,
            'superTrack' = super_name,
            'bigDataUrl' = basename(outfile),
            'shortLabel' = track_name,
            'longLabel'  = paste(track_name,' '),
            'visibility' = 'full',
            'type'       = type,
            'color'      = '34,139,34',
            'priority'   = paste(i,tn,sep='.'),
            'autoScale'  = 'on'
            )
            cat(join_string(track_description,'\t'), file=file_trackDb, append=TRUE)
        }
    
    }
    cat('done', file=file.path(path_hub,'done.txt'))
}


Make_UCSC_HUB(
  path_hub    = snakemake@params[['path_hub']],
  genome_name = snakemake@params[['genome_name']],
  hub         = snakemake@params[['hub']], 
  paths       = snakemake@params[['paths']])
