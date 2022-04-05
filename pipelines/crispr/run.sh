#!/usr/bin/bash

settings=$1
nodes=$2
additional_args=$3

export PYTHONPATH='~/.conda/envs/crispr_dart/lib/python3.5/site-packages/'
echo "ARGS: ${additional_args}"
SRCDIR="." # provide full path to this folder if running in a different directory
snakefile="${SRCDIR}/snakefile.py"
snakemake=`which snameke`
${snakemake}  ${additional_args} -p --configfile ${settings}  --snakefile ${snakefile} -j ${nodes} 
