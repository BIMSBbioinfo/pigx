

# ---------------------------------------------------------------------------- #
table_RleList = function(x)
{
    library(S4Vectors)
    nbins = max(max(x)) + 1L
    data = lapply(x,
                   function(xi)
                       S4Vectors:::tabulate2(runValue(xi) + 1L, nbins,
                                             weight=runLength(xi)))
    ans = matrix(unlist(data, use.names=FALSE),
                  nrow=length(x), byrow=TRUE)
    dimnames(ans) = list(names(x), 0:(nbins-1L))
    class(ans) = "table"
    ans
}


# ---------------------------------------------------------------------------- #
RleSumAny = function (e1, e2)
{
    library(chipseq)
    len = length(e1)
    stopifnot(len == length(e2))
    x1 = runValue(e1); s1 = cumsum(runLength(e1))
    x2 = runValue(e2); s2 = cumsum(runLength(e2))
    .Call("rle_sum_any",
          as.integer(x1), as.integer(s1),
          as.integer(x2), as.integer(s2),
          as.integer(len),
          PACKAGE = "chipseq")
}


# ---------------------------------------------------------------------------- #
Append_List_Element = function(l, name, value){
    l[[name]] = c(l[[name]], value)
    l
}

# ---------------------------------------------------------------------------- #
Summarize_Statistics_List = function(
    lout
){
    lsum = list()
    lsum$CovHistAll = NULL

    message('CovHistAll ...')
        for(l in 1:length(lout$CovHist)){
            tempCovHist = cbind(as.numeric(names(lout$CovHist[[l]])),as.numeric(lout$CovHist[[l]]))
            tempCovHist = merge(lsum$CovHistAll,tempCovHist,by.x=0,by.y=1,all=TRUE,sort=FALSE)
            CovSums    = data.frame(Depth=rowSums(tempCovHist[,-1,drop=FALSE],na.rm=TRUE),row.names=tempCovHist[,1])
            CovHistAll = CovSums
        }
        tempCovHistAll = as.numeric(CovHistAll[,1])
        names(tempCovHistAll) = rownames(CovHistAll)
        lsum$CovHistAll = tempCovHistAll

    message('ShiftMat ...')
        lsum$ShiftsAv = rowMeans(do.call(cbind, lout$ShiftMat), na.rm=TRUE)

    message('ShiftsCorAv ...')
        lsum$ShiftsCorAv = rowMeans(do.call(cbind, lout$ShiftMatCor), na.rm=TRUE)

    message('PosAny ...')
        lsum$PosAny = with(lout, unname((sum(NegAny)))+unname((sum(PosAny))))

    message('SSD ...')
        lsum$SSDAv = mean(unlist(lout$SSD))

    message('readlength ...')
        lsum$readlength = lout$readlength

    return(lsum)
}
