##' construct SummarizedExperiment
##'
##' construct SummarizedExperiment
##' the intersection of rows and columns will be used ...
##' @title construct SummarizedExperiment
##' @param assaysData data matrix
##' @param colData phenotype/metadata
##' @param rowData feature information
##' @param author optional author information
##' @param dbVersion optional version of metadatabase
##' @param note optional note will be added to exptData
##' @return object of class 'SummarizedExperiment'
##' @author mvaniterson
##' @import SummarizedExperiment
##' @importFrom S4Vectors SimpleList DataFrame
##' @importFrom utils packageVersion
##' @export
makeSE <- function(assaysData, colData, rowData, author=NULL, dbVersion=NULL,
                   note=NULL){
    ##simple constructor

    exptData <- SimpleList(creationDate = date(),
                           author = author,
                           dbVersion = dbVersion,
                           BBMRIomicsVersion = as.character(packageVersion("BBMRIomics")),
                           note = note)

    ##match rowData with assaysData
    rid <- match(rownames(assaysData), names(rowData))
    rowData <- rowData[rid]

    ##match colData with assaysData
    cid <- match(colnames(assaysData), rownames(colData))
    colData <- colData[cid,]
    colData <- droplevels(colData)

    SummarizedExperiment(assays=SimpleList(data=assaysData),
                         rowRanges=rowData,
                         colData=DataFrame(colData),
                         metadata=exptData)
}

is.unrelated <- function(cData){
    
    relations <- cData$relation_uuid
    unrelated <- TRUE | logical(nrow(cData))
    names(relations) <- names(unrelated) <- cData$uuid

    relations <- relations[!is.na(relations)]

    if(length(relations) < 1)
        return(invisible(unrelated))
    
    keep <- drop <- c()
    for(i in 1:length(relations)) {
        keeping <- names(relations[i])
        dropping <- unlist(strsplit(relations[i], ","))
        if(!(keeping %in% drop)) {
            keep <- c(keep, keeping)
            drop <- c(drop, dropping)
        }
    }

    unrelated[names(unrelated) %in% drop] <- FALSE

    return(invisible(unrelated))
}

