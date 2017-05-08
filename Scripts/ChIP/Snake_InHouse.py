
import glob
import os
import re


CHECK = os.environ['MYLIB'] + '/Perl/CheckCoord.pl'
files = [file for file in glob.glob('./**/*.bam', recursive=True)]
NAMES  = [re.sub('.bam','',os.path.basename(file)) for file in files]
MAPPED_DIR = 'Mapped/Bowtie2'
GENOME = ''
QC_PATH = ''
ex_nams = dict(zip(NAMES,NAMES))
print(ex_nams)

# ----------------------------------------------------------------------------- #
include: '/home/vfranke/Projects/CBirchmeier_Neuro/Scripts/Snake/Snake_Functions.py'
localrules: makelinks

rule all:
	input:
		expand(os.path.join(MAPPED_DIR, "{name}", "{name}.bw"),  name=NAMES) +
		expand(os.path.join(MAPPED_DIR, "Tracks",  "{ex_name}.bw"),  ex_name=ex_nams.keys())
