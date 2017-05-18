# Starting FASTQ files
export FQ1=sampleB.pe1.fq
export FQ2=sampleB.pe2.fq

# The names of the random subsets you wish to create
export FQ1SUBSET=1.rand.2k.fq
export FQ2SUBSET=2.rand.2k.fq

# How many random pairs do we want?
export N=2000

# paste the two FASTQ such that the 
# header, seqs, seps, and quals occur "next" to one another
  paste $FQ1 $FQ2 | \
# "linearize" the two mates into a single record.  Add a random number to the front of each line
  awk 'BEGIN{srand()}; {OFS="\t"; \
                        getline seqs; getline sep; getline quals; \
                        print rand(),$0,seqs,sep,quals}' | \
# sort by the random number
  sort -k1,1 | \
# grab the first N records
  head -n $N | \
# Convert the stream back to 2 separate FASTQ files.
  awk '{OFS="\n"; \
        print $2,$4,$6,$8 >> ENVIRON["FQ1SUBSET"]; \
        print $3,$5,$7,$9 >> ENVIRON["FQ2SUBSET"]}'

