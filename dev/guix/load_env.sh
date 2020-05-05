export GUIX_PROFILE=$PWD/dev/guix/.guix-profile 
source ${GUIX_PROFILE}/etc/profile
export GUIX_LOCPATH=${GUIX_PROFILE}/lib/locale
export PIGX_UNINSTALLED=1



# ======  IF YOU WANT THE STANDARDIZED PACKAGE: 
# /gnu/var/guix/profiles/custom/pigx/activate 
# it cleans your environment and puts just the latest version of `pigx` onto your PATH.  It also sets up locales and SSL certs, and indicates in the shell prompt that you are in the "pigx" environment.
# Inside of this environment please use `pigx bsseq`, `pigx rnaseq`, `pigx chipseq` or `pigx scrnaseq`.  Do not set PIGX_UNINSTALLED or anything else.
# This should just work.


