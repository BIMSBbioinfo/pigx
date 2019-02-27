# ---------------------------------------------------------------------------- #
options = commandArgs(trailingOnly=TRUE)
source(file.path(options[2],'/Argument_Parser.R'))
argv = Parse_Arguments('ChIPQC')

# ---------------------------------------------------------------------------- #
#' ChIPQC
#' Current implementation works only on the largest chromosome - it does not support multiple chromosomes'
#'
#' @param bamfile
#' @param nucleotide_frequency
#' @param outfile    - location of the output file
#' @param scriptdir  - location of R scripts
#' @param config
#' @param sample_name
#' @param use_longest_chr
#' @param shift_window
#'
#' @return saves RDS object with an annotated GRanges object
# TO DO - make the function work for no mappd reads


ChIPQC = function(
    bamfile              = NULL,
    logfile              = NULL,
    nucleotide_frequency = NULL,
    annotation           = NULL,
    outfile              = NULL,
    scriptdir            = NULL,
    library_type         = NULL,
    sample_name          = NULL,
    use_longest_chr      = NULL,
    shift_window         = 400,
    threads              = 1
){

    # ------------------------------------------------------------------------ #
    # checks for default arugments
    deflist = as.list(formals(ChIPQC))
    arglist = as.list(match.call)
    arg.ind = names(deflist) %in% names(arglist)
    if(any(arg.ind))
        stop(paste(paste(names(arglist)[arg.ind], collapse=','),'not defined'))

    # ------------------------------------------------------------------------ #
    suppressPackageStartupMessages(require(dplyr))
    suppressPackageStartupMessages(require(Rsamtools))
    suppressPackageStartupMessages(require(GenomicRanges))
    suppressPackageStartupMessages(require(GenomicAlignments))
    suppressPackageStartupMessages(require(BiocParallel))
    source(file.path(scriptdir, 'ChIPQC_Functions.R'))
    source(file.path(scriptdir, 'Functions_Helper.R'))
    register(MulticoreParam(workers = threads))
    # ------------------------------------------------------------------------ #
    annotation = readRDS(annotation)$genomic_annotation

    # check for library type - determines the reading function
    if(!library_type %in% c('single','paired'))
      stop('supported library types are single or paired')

    cnts.stat = readRDS(logfile)
    mapped.total = cnts.stat[cnts.stat$sample_name == sample_name,]
    mapped.total = mapped.total$value[mapped.total$stat == 'mapped.total']
    # compensates for double number of paired reads
    # normalizes to fragment number
    cnt_factor = ifelse(library_type == 'single', 1, 2)
    mapped.total = mapped.total/cnt_factor

    chr_lengths = scanBamHeader(bamfile)[[1]]$targets
    chr_lengths = sort(chr_lengths, decreasing=TRUE)
    use_longest_chr = TRUE
    if(use_longest_chr == TRUE || use_longest_chr == 'TRUE')
        chr_lengths = head(chr_lengths, 1)

    # removes chromosomes which are shorter than shift_size
    if(length(chr_lengths[chr_lengths < shift_window]) > 0){
        message("Removing ",length(chr_lengths[chr_lengths < shift_window]),
                "chromosomes with length less than cross-coverage shift")
        chr_lengths = chr_lengths[!chr_lengths < shift_window]
        if(length(chr_lengths) == 0)
          stop('shift_size was larger than all chromoromes')
    }

    # checks whether the bamfile is indexed
    if(file.exists(bamfile) & length(index(BamFile(bamfile))) == 0){
        message("Creating index for ",bamfile)
          indexBam(bamfile)
        message("..done")
    }


    # list for collecting summary stats
    lout = list(
        ShiftMat      = NULL,
        ShiftMatCor   = NULL,
        CovHist       = NULL,
        Cov           = NULL,
        SSD           = NULL,
        PosAny        = NULL,
        NegAny        = NULL,
        SSDBL         = NULL,
        readlength    = NULL
    )

    # loops through the chromosomes and collects summary statistics
    for(k in 1:length(chr_lengths)){
#
        chrname = names(chr_lengths)[k]
        Param = ScanBamParam(
                    which=GRanges(seqnames=chrname,IRanges(start=1,end=unname(chr_lengths[chrname])-shift_window)))

        # read the reads
        if(library_type == 'single')
            temp  = granges(readGAlignments(bamfile, param=Param), use.mcols = TRUE)

        if(library_type == 'paired')
            temp  = unlist(range(grglist(readGAlignmentPairs(bamfile, param=Param))))

        if(k == 1){
            lout$readlength= width(temp)[1:(min(length(temp),1000))]
        }

        granges    = granges(temp, use.mcols=TRUE, use.names=TRUE)
        Sample_GIT = GNCList(granges)

        # calculates the coverage
        lout$Cov   = coverage(Sample_GIT,width=unname(chr_lengths[k]))

        # coverage distribution
        lout = Append_List_Element(lout, 'CovHist', list(colSums(table_RleList(lout$Cov))))

        # coverage standard deviation
        lout = Append_List_Element(lout, 'SSD', sd(lout$Cov)[chrname])

        PosCoverage = coverage(Sample_GIT[strand(Sample_GIT)=="+"],width=unname(chr_lengths[k]))[[chrname]]
        NegCoverage = coverage(Sample_GIT[strand(Sample_GIT)=="-"],width=unname(chr_lengths[k]))[[chrname]]

        lout = Append_List_Element(lout, 'PosAny', sum(runLength((PosCoverage[PosCoverage > 1]))))
        lout = Append_List_Element(lout, 'NegAny', sum(runLength((NegCoverage[NegCoverage > 1]))))

        # Cross correlation
        ShiftsTemp = shiftApply(seq(1, shift_window, 10),PosCoverage,NegCoverage,RleSumAny)
        lout = Append_List_Element(lout, 'ShiftMat', setNames(list(ShiftsTemp), chrname))

        ShiftsCorTemp = shiftApply(seq(1,shift_window, 10),PosCoverage,NegCoverage,cor)
        lout = Append_List_Element(lout, 'ShiftMatCor', setNames(list(ShiftsCorTemp), chrname))


        # Read duplication - library complexity calculation
        granges_unique = unique(granges)
        lout = Append_List_Element(lout, 'duplicated_reads',  length(granges_unique))

        
        cov_uniq = coverage(granges_unique, width=unname(chr_lengths[k]))
        cov_ind = cov_uniq > 0
        lout = Append_List_Element(lout, 'libComplexity', round(cov_uniq[cov_ind]/ lout$Cov[cov_ind],2))

        lout$annot = data.frame(annot = AnnotateRanges(temp, annotation)) %>%
              group_by(annot) %>%
              summarize(N = n()) %>%
              mutate(total = sum(N)) %>%
              mutate(freq  = round(N/total,3)) %>%
              dplyr::select(annot, freq) %>%
              mutate(sample_name = sample_name)

    }

    lsum = Summarize_Statistics_List(lout)

    tilling_windows = readRDS(nucleotide_frequency)
    cnts = summarizeOverlaps(tilling_windows, bamfile)
    cnts = DataFrame(assays(cnts)[[1]])
    cnts$norm = round(cnts[[1]]*(1e6/mapped.total),2)
    tilling_windows$samplecounts = cnts
    lsum$tilling_windows = tilling_windows

    saveRDS(lsum, file = outfile)
}



# ---------------------------------------------------------------------------- #
ChIPQC(
  bamfile              = argv$input[['bamfile']],
  logfile              = argv$input[['logfile']],
  nucleotide_frequency = argv$input[['nucleotide_frequency']],
  annotation           = argv$input[['annotation']],
  outfile              = argv$output[['outfile']],
  scriptdir            = argv$params[['scriptdir']],
  library_type         = argv$params[['library_type']],
  shift_window         = argv$params[['shift_window']],
  sample_name          = argv$params[['sample_name']],
  use_longest_chr      = argv$params[['use_longest_chr']],
  threads              = argv$params[['threads']]
)
