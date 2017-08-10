#!/usr/bin/env python3.5
import os, sys, json



# -------------------------------------------------------------------------------
# Without the main sentinel, the code would be executed even if the script were imported as a module.

def main(argv):

    """
    The main function here passes all of the parameters to a cluster configuration file in JSON format 
    it is called from the main shell script once for each rule that requires parameter specification
    (plus once more for the default configuration.) when "comma"==0, we should be at the last item in the list.
    """


    if len(argv)!= 7:
      print("ERROR: "+str(len(argv))+" command line arguments supplied to create_cluster. 6 are required.")
      exit()
    rulename  = argv[1]
    nthreads  = argv[2]
    q         = argv[3]
    MEM       = argv[4]
    h_stack   = argv[5]
    comma     = int(argv[6])

 
    print("    \""+rulename+"\":" )
    print("    {" )
    print("        \"nthreads\":"+nthreads+"," )
    print("        \"MEM\":\""+MEM+"\"," )
    print("        \"q\":\""+q+"\"," )
    print("        \"h_stack\":\""+h_stack+"\"" )

    if ( comma == 1 ) :
      print("    },")
    else:
      print("    }")
      print("}")
  
  
if __name__ == "__main__":
    main(sys.argv)

   
   
   
   
   
   
  
   
   
   
   
   
   
   
  
   
   
   
   
   
   
   
   
   
   
