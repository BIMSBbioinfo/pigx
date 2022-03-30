# PiGx RNAseq Pipeline.
#
# Copyright © 2022 Bora Uyar <bora.uyar@mdc-berlin.de>
#
# This file is part of the PiGx RNAseq Pipeline.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Check the contents of the annotation files (e.g. GTF/Fasta)
# for potential errors to abort early. 

args <- commandArgs(trailingOnly = TRUE)

gtfFile <- args[1] # path to gene annotations in GTF format
cDNAfastaFile <- args[2] # path to transcriptome sequences in fasta format 
genomeFastaFile <- args[3] # path to DNA sequences in fasta format
outDir <- args[4] # where to save the annotation stats

message(date(), " Checking annotation files for potential issues")

# Check to see if the input GTF file is parseable.
message(date()," => Checking to see if GTF file can be properly imported from here: ",gtfFile)
gtf <- rtracklayer::import.gff(gtfFile)
message(date()," => Imported ",length(gtf)," features from the GTF file")
if(length(gtf) == 0) {
  stop("ERROR: GTF file seems to be empty here:",gtfFile)
} 

# check if gene_id field is accessible
if(!'gene_id' %in% colnames(GenomicRanges::mcols(gtf))){
  stop("ERROR: Can't access gene_id column in the GTF File, check the format of your GTF file")
}
if(!'transcript_id' %in% colnames(GenomicRanges::mcols(gtf))){
  stop("ERROR: Can't access transcript_id column in the GTF File, check the format of your GTF file")
}

message(date(), " => GTF file contains annotations for ",length(unique(gtf$gene_id)), " genes ",
        "and ",length(unique(gtf$transcript_id)), " transcripts")

# Check to see if the cDNA fasta file can be imported 
cDNA <- Biostrings::readDNAStringSet(cDNAfastaFile, format = 'fasta')
if(length(cDNA) == 0){
  stop("ERROR: cDNA fasta file seems to be empty")
}
message(date(), " => Imported ",length(cDNA), " transcripts from the cDNA fasta file at", cDNAfastaFile)
# Check to see if the transcript ids in the GTF file match the transcript ids in the cDNA Fasta file
# if not, give a warning how this will affect id matching for Salmon and give a suggestion on 
# how to update the cDNA file 

matched <- length(intersect(names(cDNA), unique(gtf$transcript_id)))
message(date(), " => Number of transcript ids matching between the GTF file and the cDNA fasta file: ",matched)
if(matched < 0.1 * length(unique(gtf$transcript_id))) {
  warning(date(), " => Couldn't match the transcript ids between the GTF file and the cDNA fasta file. 
          This will cause problems in getting gene-level expression estimates when running Salmon.  
          However, the transcript-level expression estimation won't suffer. 
          Possible reason is the source databases of these annotations are different. 
          Another possible reason is the additional annotations on the transcript IDs in the 
          fasta file: example entry for the transcript 'ENST00000202017' in GTF file is 
          'ENST00000202017.5 cdna chromosome:GRCh38:20:31944342:31952092 ...' in the Fasta File. 
          If using annotations from ENSEMBL, the fasta file can be cleaned up with the following 
          `sed` command: ",
          "sed 's/\\(ENST[0-9]*\\)\\.[0-9]*.*/\\1/g' <cDNA_file_path> > cDNA_file.cleaned.fasta")
}


# Check to see if the genome fasta can be imported and 
# see if the chromosome names in the genome fasta are on 
# the same convention as used in the GTF file
DNA <- Biostrings::readDNAStringSet(genomeFastaFile, format = 'fasta')
message(date(), " => Imported ",length(DNA)," chromosomes/contigs from the DNA fasta files at",genomeFastaFile)

# abort if the chromosome names don't match in DNA and GTF 
match <- intersect(GenomeInfoDb::seqlevelsStyle(DNA), GenomeInfoDb::seqlevelsStyle(gtf))
message(date(), " Chromosome naming convention for the genome:",
        paste(GenomeInfoDb::seqlevelsStyle(DNA), collapse = ' '))
message(date(), " Chromosome naming convention for the GTF file:", 
        paste(GenomeInfoDb::seqlevelsStyle(gtf), collapse = ' '))
if(length(match) == 0) {
  stop("ERROR: the chromosome naming conventions in the DNA file and the GTF file don't agree.\n",
       "Example names from DNA file:",paste(head(unique(names(DNA))), collapse = ' '),
       "\nExample names from GTF file:",paste(head(unique(GenomicRanges::seqnames(gtf))), collapse = ' '),
       "\nUsing genome and annotation files from the same source (e.g. the ENSEMBL database)", 
       " is highly recommended.")
}

# print out some stats about the input annotations

input_stats <- rbind(data.frame('Type' = 'GTF', 
                                'Description' = 'Number of genes',
                                'Value' = length(unique(gtf$gene_id))),
                     data.frame('Type' = 'GTF', 
                                'Description' = 'Number of transcripts', 
                                'Value' = length(unique(gtf$transcript_id))),
                     data.frame('Type' = 'GTF', 
                                'Description' = 'Naming convention', 
                                'Value' = paste(GenomeInfoDb::seqlevelsStyle(gtf), collapse = ' ')),
                     data.frame('Type' = 'DNA', 
                                'Description' = 'Number of chromosomes', 
                                'Value' = length(DNA)), 
                     data.frame('Type' = 'DNA', 
                                'Description' = 'Naming convention', 
                                'Value' = paste(GenomeInfoDb::seqlevelsStyle(DNA), collapse = ' ')), 
                     data.frame('Type' = 'cDNA', 
                                'Description' = 'Number of transcripts', 
                                'Value' = length(cDNA)), 
                     data.frame('Type' = 'cDNA', 
                                'Description' = 'Mean length of transcripts', 
                                'Value' = mean(lengths(cDNA))), 
                     data.frame('Type' = 'cDNA vs GTF', 
                                'Description' = 'Number of matching transcript identifiers', 
                                'Value' = length(intersect(names(cDNA), unique(gtf$transcript_id)))))

write.table(input_stats, file = file.path(outDir, 'input_annotation_stats.tsv'), quote = F, row.names = F)

message(date(), " => Finished checking annotation files.")




