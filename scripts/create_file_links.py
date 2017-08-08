#!/usr/bin/env python3.5

import os, sys, json
from helpers import *

########################################################

# Without the main sentinel, the code would be executed even if the script were imported as a module.
def main(argv):
    """
    The main function to create file links from a table sheet
    (and optionally paths to tools in the JSON file).
    """

    if len(argv) != 2:
      print("Error in create_file_links.py: script requires one command line argument specifying the config file. Exiting.")
      exit()

    configfile = argv[1]   
    if not os.path.isfile(configfile):
      print("ERROR: create_file_links is looking for a config file that doesn't exist")
    else:
      json_datastream=open(configfile)
      config=json.load(json_datastream)

      path_SOURCE  = config['PATHIN'];  
      if not path_SOURCE.endswith("/"): 
        path_SOURCE = path_SOURCE+"/"
      
      path_OUT     = config['PATHOUT']; 
      if not path_OUT.endswith("/"):    
        path_OUT    = path_OUT+"/"

      for s in config['SAMPLES']:
        flist = config['SAMPLES'][s]['files']

        # -------- single-end case --------------
        if len(flist) == 1: # Single-end

          if not flist[0].endswith(".gz"):
            print("error: input files must be gzipped. this job will fail. todo: subsequent versions will handle unzipped .fq or .bz2.")

          linkname = config['SAMPLES'][s]['SampleID']+".fq.gz"

          #--- now make the single link to this one file: 
          makelink(path_SOURCE+flist[0], path_OUT+"path_links/input/", linkname )

        # -------- paired-end case --------------
        elif len(flist) == 2:

          if ( not flist[0].endswith(".gz") ) or ( not flist[1].endswith(".gz") ) :
            print("error: input files must be gzipped. this job will fail. todo: subsequent versions will handle unzipped .fq or .bz2.")
        
          linkname_1 = config['SAMPLES'][s]['SampleID']+"_1.fq.gz"
          linkname_2 = config['SAMPLES'][s]['SampleID']+"_2.fq.gz"
          
          makelink(path_SOURCE+flist[0], path_OUT+"/path_links/input/" , linkname_1)
          makelink(path_SOURCE+flist[1], path_OUT+"/path_links/input/" , linkname_2)

      
if __name__ == "__main__":
    main(sys.argv)  
 

