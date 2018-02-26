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
# bamfile  = '/data/local/Projects/pigx/pigx_chipseq/Data/Analyzed/Mapped/Bowtie/mm_tet1_d1_h3k9me3_br2/mm_tet1_d1_h3k9me3_br2.sorted.bam'
# nucleotide_frequency = '/data/local/Projects/pigx/pigx_chipseq/Data/Analyzed/Bowtie2_Index/mm9/mm9.NucleotideFrequency.GRanges.rds'
# sample_sheet = yaml::yaml.load_file('/home/vfranke/Projects/AAkalin_pigx/sample_sheet.yaml')
# sample_name = 'mm_tet1_d1_h3k9me3_br2'
# scriptdir = 'scripts'



ChIPQC = function(
    bamfile              = NULL,
    logfile              = NULL,
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
        #     colnames(ShiftMattemp) = names(chr_lengths)[k]
        #     lout[['ShiftMat']] = cbind(lout[['ShiftMat']],ShiftMattemp)
        #
        # }else{
            if(k == 1){
                lout$readlength= width(temp)[1:(min(length(temp),1000))]
            }
            Sample_GIT = GNCList(granges(temp, use.mcols=TRUE, use.names=TRUE))
            lout$Cov   = coverage(Sample_GIT,width=unname(chr_lengths[k]))

            lout = Append_List_Element(lout, 'CovHist', list(colSums(table_RleList(lout$Cov))))
            lout = Append_List_Element(lout, 'SSD', sd(lout$Cov)[chrname])

            PosCoverage = coverage(Sample_GIT[strand(Sample_GIT)=="+"],width=unname(chr_lengths[k]))[[chrname]]
            NegCoverage = coverage(Sample_GIT[strand(Sample_GIT)=="-"],width=unname(chr_lengths[k]))[[chrname]]

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
  outfile              = argv$output[['outfile']],
  scriptdir            = argv$params[['scriptdir']],
  sample_sheet         = argv$params[['sample_sheet']],
  sample_name          = argv$params[['sample_name']],
  use_longest_chr      = argv$params[['use_longest_chr']]
)
