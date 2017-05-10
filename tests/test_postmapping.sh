#!/bin/bash
## This small test is used to check wether the post-mapping steps work
## 
## If there are any R related issues, check wether the dependencies listed in 
## rules/postmapping.rules are installed.
##
## please run from the same directory where BSseq_pipeline.py is located

echo "Testing post mapping steps"

snakemake --snakefile BSseq_pipeline.py --configfile new.config.json annotation/sampleA.pe1_trimmed_bismark_bt2.deduplicated.sorted_ce10_annotation.nb.html  --forceall -n