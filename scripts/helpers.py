import os

def makelink( actual_file, link_location, linkname):
  if not os.path.isfile(actual_file):
    print("ERROR: linking to non-existent file: "+actual_file ); status = -1
  elif not os.path.isdir(link_location):
    print("ERROR: PATHOUT and/or subdirectory does not exist for linking: "+link_location); status = -1
  else:
    status = os.system( "ln -sn "+actual_file+"  "+link_location+linkname+" 2>/dev/null" )
    #-- create a link only if non previously exists. suppress warnings (silent) if it's already there.
  return status
