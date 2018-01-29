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

# basepath = './'
# bamfile  = file.path(basepath, 'Tests/out/Mapped/Bowtie/ChIP1/ChIP1.sorted.bam')
# outfile  = file.path(basepath, 'Tests/out/Analysis/RDS/')
# scriptdir   = file.path(basepath, 'scripts')
# genomic_windows = file.path(basepath)
# nucleotide_frequency = file.path(basepath,'Tests/out/Bowtie2_Index/hg19/hg19.NucleotideFrequency.GRanges.rds')


ChIPQC = function(
    bamfile              = NULL,
    nucleotide_frequency = NULL,
    outfile              = NULL,
    scriptdir            = NULL,
    sample_sheet         = NULL,
    sample_name          = NULL,
    use_longest_chr      = NULL,
    shift_window         = 400
){

    # ------------------------------------------------------------------------ #
    # checks for default arugments
    deflist = as.list(formals(ChIPQC))
    arglist = as.list(match.call)
    arg.ind = names(deflist) %in% names(arglist)
    if(any(arg.ind))
        stop(paste(paste(names(arglist)[arg.ind], collapse=','),'not defined'))

    # ------------------------------------------------------------------------ #
    suppressPackageStartupMessages(require(Rsamtools))
    suppressPackageStartupMessages(require(GenomicRanges))
    suppressPackageStartupMessages(require(GenomicAlignments))
    source(file.path(scriptdir, 'ChIPQC_Functions.R'))
    # ------------------------------------------------------------------------ #


    # check for library type - determines the reading function
    library_type = sample_sheet$samples[[sample_name]]$library
    if(!library_type %in% c('single','paired'))
      stop('supported library types are single or paired')

    ChrLengths = scanBamHeader(bamfile)[[1]]$targets
    ChrLengths = sort(ChrLengths, decreasing=TRUE)
    use_longest_chr = TRUE
    if(use_longest_chr == TRUE || use_longest_chr == 'TRUE')
        ChrLengths = head(ChrLengths, 1)

    # removes chromosomes which are shorter than shift_size
    if(length(ChrLengths[ChrLengths < shift_window]) > 0){
        message("Removing ",length(ChrLengths[ChrLengths < shift_window]),
                "chromosomes with length less than cross-coverage shift")
        ChrLengths = ChrLengths[!ChrLengths < shift_window]
        if(length(ChrLengths) == 0)
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
    for(k in 1:length(ChrLengths)){
#
        chrname = names(ChrLengths)[k]
        Param = ScanBamParam(
                    which=GRanges(seqnames=chrname,IRanges(start=1,end=unname(ChrLengths[chrname])-shift_window)))

        if(library_type == 'single')
            temp  = granges(readGAlignments(bamfile, param=Param), use.mcols = TRUE)

        if(library_type == 'paired')
            temp  = unlist(range(grglist(readGAlignmentPairs(bamfile, param=Param))))
        #
        # if(length(temp) == 0){
        #
        #     setnames = c('Chr_SSD','SSDBL','PosAny','NegAny')
        #     for(i in setnames){
        #         v = 0; names(v) =  chrname
        #         lout[[i]] = c(lout[[i]], v)
        #     }
        #     emptyChr_SSD = 0
        #
        #     ShiftMattemp = matrix(rep(0,(shiftWindowEnd-shiftWindowStart)+1),ncol=1)
        #     colnames(ShiftMattemp) = names(ChrLengths)[k]
        #     lout[['ShiftMat']] = cbind(lout[['ShiftMat']],ShiftMattemp)
        #
        # }else{
            if(k == 1){
                lout$readlength=round(median(width(temp)))
                lout$readlength_dist = hist(width(temp), breaks=50)
            }
            Sample_GIT = GNCList(granges(temp, use.mcols=TRUE, use.names=TRUE))
            lout$Cov   = coverage(Sample_GIT,width=unname(ChrLengths[k]))

            lout = Append_List_Element(lout, 'CovHist', list(colSums(table_RleList(lout$Cov))))
            lout = Append_List_Element(lout, 'SSD', sd(lout$Cov)[chrname])

            PosCoverage = coverage(Sample_GIT[strand(Sample_GIT)=="+"],width=ChrLengths[k])[[chrname]]
            NegCoverage = coverage(Sample_GIT[strand(Sample_GIT)=="-"],width=ChrLengths[k])[[chrname]]

            lout = Append_List_Element(lout, 'PosAny', sum(runLength((PosCoverage[PosCoverage > 1]))))
            lout = Append_List_Element(lout, 'NegAny', sum(runLength((NegCoverage[NegCoverage > 1]))))

            ShiftsTemp = shiftApply(seq(shift_window),PosCoverage,NegCoverage,RleSumAny)
            lout = Append_List_Element(lout, 'ShiftMat', setNames(list(ShiftsTemp), chrname))

            ShiftsCorTemp = shiftApply(seq(shift_window),PosCoverage,NegCoverage,cor)
            lout = Append_List_Element(lout, 'ShiftMatCor', setNames(list(ShiftsCorTemp), chrname))

            # }
    }

    lsum = Summarize_Statistics_List(lout)


    tilling_windows = readRDS(nucleotide_frequency)
    cnts = summarizeOverlaps(tilling_windows, bamfile)
    tilling_windows$samplecounts = DataFrame(assays(cnts)[[1]])
    lsum$tilling_windows = tilling_windows

    saveRDS(lsum ,file = outfile)
}



# ---------------------------------------------------------------------------- #
ChIPQC(
  bamfile              = argv$input[['bamfile']],
  nucleotide_frequency = argv$input[['nucleotide_frequency']],
  outfile              = argv$output[['outfile']],
  scriptdir            = argv$params[['scriptdir']],
  sample_sheet         = argv$params[['sample_sheet']],
  sample_name          = argv$params[['sample_name']],
  use_longest_chr      = argv$params[['use_longest_chr']]
)
