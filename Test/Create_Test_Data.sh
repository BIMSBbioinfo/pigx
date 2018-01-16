inpath='/data/local/vfranke/AAkalin_PIX/scRNA/001-42802245'

genome='/data/local/vfranke/AAkalin_PIX/scRNA/Annotation/hg19_mm9/STAR_INDEX'
infile=(`find $inpath | grep fastq | grep -v Test | grep R2`)
outfile=$inpath/Test
STAR='/home/vfranke/bin/Software/Mapping/STAR-2.5.3a/bin/Linux_x86_64_static/STAR'

$STAR --genomeDir $genome --readFilesIn ${infiles[@]} --runThreadN 16 --genomeLoad LoadAndKeep --outFilterMultimapNmax 10 --outFileNamePrefix $outfile --outReadsUnmapped Fastx --outFilterMismatchNoverLmax 0.05 --readFilesCommand 'zcat'