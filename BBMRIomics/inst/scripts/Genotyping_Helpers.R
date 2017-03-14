.rectangular <- function(x, y, verbose){

    message("Running `rectangular` IBS algorithm!")

    if((nrow(x) != nrow(y)))
        stop("Dimension mismatch!")

    na.rm <- sum(is.na(x)) > 0 | sum(is.na(y)) > 0

    ##calculates mean and variance of IBS between all pairs of x and y
    N <- ncol(x)
    M <- ncol(y)
    K <- nrow(x)
    indices <- 1:N
    mn <- s2 <- numeric(length=N*M)
    for(j in 1:M) {
        ibs <- 2 - abs(x - y[,j])
        mn[indices + N*(j-1)] <- colMeans(ibs, na.rm=na.rm)
        s2[indices + N*(j-1)] <- colVars(ibs, na.rm=na.rm)
        ##slower
        ##tabs <- apply(ibs, 2, tabulate)
        ##mn[indices + N*(j-1)] <- crossprod(c(1,2), tabs)/K
        ##xs2[indices + N*(j-1)] <- (crossprod(c(1,4), tabs) - K*mn[indices + N*(j-1)]^2)/(K-1)
        if(verbose & (j %% floor(M/10) == 0 | j == 0))
            message(j*N, " of ", N*M, " (", round(100*j/M, 2), "%) ...")
    }
    data.frame(mean=mn,
               var=s2,
               colnames.x=rep(colnames(x), M),
               colnames.y= rep(colnames(y), each=N))
}

.square <- function(x, y, verbose){

    message("Running `square` IBS algorithm!")

    if((ncol(x) != ncol(y)) | (nrow(x) != nrow(y)))
        stop("Dimension mismatch!")

    na.rm <- sum(is.na(x)) > 0 | sum(is.na(y)) > 0

    ##calculates mean and variance of IBS between all unique pairs of x (and y)
    N <- ncol(x)
    K <- nrow(x)
    k <- 0
    mn <- s2 <- numeric(length=N*(N+1)/2)
    for(j in 1:N) {
        indices <- j:N
        k <- k + length(indices)
        ibs <- 2 - abs(x[, indices, drop=FALSE] - y[,j])
        mn[k-length(indices) + indices-(j-1)] <- colMeans(ibs, na.rm=na.rm)
        s2[k-length(indices) + indices-(j-1)] <- colVars(ibs, na.rm=na.rm)
        ##slower
        ##tabs <- apply(ibs, 2, tabulate)
        ##mn[k-length(indices) + indices-(j-1)] <- crossprod(c(1,2), tabs)/K
        ##s2[k-length(indices) + indices-(j-1)] <- (crossprod(c(1,4), tabs) - K*mn[k-length(indices) + indices-(j-1)]^2)/(K-1)

        if(verbose & (j %% 100 == 0 | j == 1))
            message(k, " of ", N*(N+1)/2, " (", round(100*k/(N*(N+1)/2), 2), "%) ...")
    }

    data.frame(mean = mn,
               var = s2,
               colnames.x = unlist(sapply(1:N, function(k) colnames(x)[k:N])),
               colnames.y = rep(colnames(y), N:1))
}


.phasing <- function(x, y, rHash) {
    ##relabel those in x according to those in y
    ##relabelling is based on the idea that snp's in x should be positively correlated with those in y
    ##robust against NA's and outliers

    ##FIX level `identical`
    identical <- names(which(unlist(mget(ls(rHash), rHash, mode = "integer", ifnotfound=list(0L)), use.names=TRUE) == 2))

    colnames.x <- gsub(":.*$", "", identical)
    colnames.y <- gsub("^.*:", "", identical)
    colnames.x <- colnames.x[!grepl("NA", colnames.y)]
    colnames.y <- colnames.y[!grepl("NA", colnames.y)]

    keep <- colnames.x %in% colnames(x) & colnames.y %in% colnames(y)
    if(sum(keep) == 0)
        stop("No overlapping samples!")

    colnames.x <- colnames.x[keep]
    colnames.y <- colnames.y[keep]

    midx <- match(colnames.x, colnames(x))
    midy <- match(colnames.y, colnames(y))

    signs <- sign(unlist(lapply(1:nrow(x), function(i) cov(x[i,midx], y[i,midy], use="complete.obs", method="spearman"))))
    for(i in 1:nrow(x)) {
        if(signs[i] < 0) {
            xi <- x[i,]
            x[i,xi==1] <- 3
            x[i,xi==3] <- 1
        }
    }
    return(x)

    ##better/faster?
    ## for(i in 1:nrow(x)) {
    ##     org <- sum(diag(table(x[i, midx], y[i, midy])))
    ##     x[i, midx] <- x[i, midx] + y[i, midy]
    ##     y[i, midy] <- x[i, midx] - y[i, midy]
    ##     x[i, midx] <- x[i, midx] - y[i, midy]
    ##     if(sum(diag(table(x[i, midx], y[i, midy]))) < org)  {
    ##         x[i, midx] <- x[i, midx] + y[i, midy]
    ##         y[i, midy] <- x[i, midx] - y[i, midy]
    ##         x[i, midx] <- x[i, midx] - y[i, midy]
    ##     }
    ## }
    ## return(x)
}


##' allele sharing based on ibs
##'
##' calculate mean variance between to vectors/matrices genotypes
##' coded as 1,2,3
##' @title allele sharing based on ibs
##' @param x genotype vector or matrix
##' @param y genotype vector or matrix
##' @param relationHash
##' @param phasing FALSE
##' @param verbose show progress
##' @return data.frame with mean and variance ibs between all pairs
##' @author mvaniterson
##' @importFrom matrixStats colVars
##' @export
alleleSharing <- function(x, y=NULL, rHash, phasing=FALSE, verbose=TRUE) {

    suppressPackageStartupMessages({
        require(matrixStats)
    })

    if(is.null(y)) {
        message("Using ", nrow(x), " polymophic SNPs to determine allele sharing.")
        data <- .square(x, x, verbose)
    } else {
        rows <- intersect(rownames(x), rownames(y))
        rId <- match(rows, rownames(x))
        x <- x[rId,]
        rId <- match(rows, rownames(y))
        y <- y[rId,]
        if(phasing)
            x <- .phasing(x, y, rHash)
        message("Using ", nrow(x), " polymophic SNPs to determine allele sharing.")
        data <- .rectangular(x, y, verbose)
    }

    pairs <- paste(data$colnames.x, data$colnames.y, sep=":")
    mId <- pairs %in% ls(rHash)
    data$relation <- "unrelated"
    data$relation[mId] <- unlist(mget(pairs[mId], rHash, mode = "character", ifnotfound=list("unrelated")), use.names=FALSE)
    data$relation <- factor(data$relation)
    invisible(data)
}

predict <- function(data, n=100, plot.it=TRUE){

    data <- droplevels(data)

    ##LDA
    suppressPackageStartupMessages({
        require(MASS)
    })

    model <- lda(relation~mean+var, data=data)

    predicted <- MASS:::predict.lda(model, data)

    data$predicted <- predicted$class

    print(table(`Predicted relation`=data$predicted, `Assumed relation`=data$relation), zero.print = ".")

    id <- which(data$predicted != data$relation)
    if(plot.it)
        plot(data[, c("mean", "var")], pch=".", cex=3, col=as.integer(data$relation), xlab="mean (IBS)", ylab="variance (IBS)")

    xp <- seq(min(data$mean), max(data$mean), length = n)
    yp <- seq(min(data$var), max(data$var), length = n)
    grid <- expand.grid(mean = xp, var = yp)
    predicted <- MASS:::predict.lda(model, grid)
    posterior <- predicted$posterior
    if(ncol(posterior) > 2) {
        for(k in 1:ncol(posterior)) {
            zp <- posterior[, k] - apply(posterior[,-k], 1, max)
            contour(xp, yp, matrix(zp, n), add = T, levels = 0, drawlabels=FALSE, lty=1, lwd=2, col="grey")
        }
    } else {
        zp <- posterior[, 2] - pmax(posterior[, 1])
        contour(xp, yp, matrix(zp, n), add = T, levels = 0, drawlabels=FALSE, lty=1, lwd=2, col="grey")
    }

    points(data[id, c("mean", "var")], pch=".", cex=3, col=as.integer(data$relation[id]))
    legend("topright", paste("assumed", levels(data$relation)), col=1:nlevels(data$relation), pch=15, bty="n")

    invisible(data[id,])
}


hashRelations <- function(relations, idx.col="id.x", idy.col="id.y", rel.col="relation_type"){
    hash <- new.env()
    keys <- paste(relations[, idx.col], relations[, idy.col], sep=":")
    values <- relations[, rel.col]

    ##are there inconsistent relations
    if(any(tapply(values, keys, function(x) length(unique(x)) > 1)))
        stop("Nonconsistent relations!")

    ##redundant relations are automatically remove!
    tmp <- mapply(assign, keys, values, MoreArgs=list(envir = hash))
    hash
}

##' convert betas to genotypes
##'
##' based on idea's from Leonard Schalkwyk (wateRmelon)
##' otherwise heterozygous hemimethylated
##' @title converts beta-values to genotypes (1, 2 and 3)
##' @param na.rm TRUE drop cp for which no clustering was observed
##' @param minSep minimal separation between clusters
##' @param minSize size of smallest cluster (in %)
##' @param betas beta matrix of probes containing SNPs
##' @return matrix with genotypes
##' @author mvaniterson
##' @importFrom stats kmeans
##' @export
beta2genotype <- function(betas, na.rm=TRUE, minSep = 0.25, minSize = 5, centers = c(0.2, 0.5, 0.8)){

    genotypes <- apply(betas, 1, .beta2genotype, minSep = minSep, minSize = minSize, centers = centers)

    genotypes <- t(genotypes)
    colnames(genotypes) <- colnames(betas)
    rownames(genotypes) <- rownames(betas)

    if(na.rm) {

        nas <- apply(betas, 2, function(x) sum(is.na(x))/length(x))

        nas <- apply(genotypes, 1, anyNA)
        genotypes <- genotypes[!nas, ]
    }

    genotypes
}

.beta2genotype <- function(x, minSep = 0.25, minSize = 5, centers = c(0.2, 0.5, 0.8)) {

    ##first check three clusters
    km <- try(kmeans(as.numeric(x), centers), silent=TRUE)
    if(!inherits(km, 'try-error')) {
        if (all(abs(rep(km$centers, 3) - rep(km$centers, each=3))[-c(1,5,9)] > minSep)) {
            if(100*min(as.numeric(table(km$cluster)))/length(x) > minSize)
                return(km$cluster)
        }
    }

    ##failed so check two clusters
    ## km <- try(kmeans(as.numeric(x), centers[-2]), silent=TRUE)
    ## if(!inherits(km, 'try-error')) {
    ##     if(all(abs(rep(km$centers, 2) - rep(km$centers, each=2))[-c(1,4)] > minSep)) {

    ##         if (km$centers[1] < centers[1] & km$centers[2] > centers[3]) ##relabelling of clusters is necessary
    ##             km$clusters[km$clusters == 2] <- 3
    ##         else if (km$centers[1] > centers[1] & km$centers[2] > centers[3])
    ##             km$clusters <- km$clusters + 1

    ##         if(100*min(as.numeric(table(km$cluster)))/length(x) > minSize)
    ##             return(km$cluster)
    ##     }
    ## }

    ##no clusters detected
    return(rep(NA, length(x)))
}


##some BBMRIomics specific helper functions
DNAmCalls <- function(cohort, DNAmFile, verbose=FALSE, maxbatch=500){

    suppressPackageStartupMessages({
        require(BBMRIomics)
        require(minfi)
        require(DNAmArray)
        require(BiocParallel)
        require(FDb.InfiniumMethylation.hg19)
    })

    samplesheets <- getView("methylationSamplesheet", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB, verbose=verbose)

    ##get location idat-files on VM
    path450k <- file.path(VM_BASE_DATA, "IlluminaHumanMethylation450k")
    samplesheets$biobank_id <- gsub("-.*$", "", samplesheets$ids)
    samplesheets$Basename <- with(samplesheets, file.path(path450k, "raw", Sentrix_Barcode,
                                                          paste(Sentrix_Barcode, Sentrix_Position, sep = "_")))

    samplesheets <- samplesheets[!duplicated(samplesheets$run_id),]

    runs <- getView("getMethylationRuns", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB, verbose=verbose)
    drop <- runs$run_id[runs$qc == "bad quality"]
    targets <- samplesheets[!(samplesheets %in% drop),]

    if(cohort == "NTR") {
        targets <- targets[targets$run_id != "9340996015_R06C02", ] ##all NA's?
    }

    if(cohort != "ALL")
        targets <- targets[targets$biobank_id == cohort,]

    load(DNAmFile)
    cpgs <- names(DNAmSNPs)

    register(MulticoreParam(6))
    if(nrow(targets) > maxbatch) {
        betas <- lapply(split(targets, 1+(1:nrow(targets))%/%maxbatch), function(targetsbatch) {
            RGset <- read.metharray.exp.par(targetsbatch, verbose=verbose)
            beta <- getBeta(RGset)
            betas <- rbind(beta[rownames(beta) %in% cpgs, ], getSnpBeta(RGset))
        })
        betas <- do.call('cbind', betas)
    } else {
        RGset <- read.metharray.exp.par(targets, verbose=verbose)
        betas <- getBeta(RGset)
        betas <- rbind(betas[rownames(betas) %in% cpgs,], getSnpBeta(RGset))
    }

    beta2genotype(betas)
}


RNACalls <- function(vcfFile, mafThr=0.025, verbose){

    suppressPackageStartupMessages({
        require(VariantAnnotation)
        require(Rsamtools)
    })

    vcf <- readVcf(vcfFile, "hg19")
    rnaCalls <- genotypeToSnpMatrix(vcf)$genotypes
    rnaCalls <- t(matrix(as.numeric(rnaCalls), nrow=nrow(rnaCalls), ncol=ncol(rnaCalls), dimnames=dimnames(rnaCalls)))

    rownames(rnaCalls) <- gsub("_.*$", "", rownames(rnaCalls))
    colnames(rnaCalls) <- gsub("\\.variant.*$", "", colnames(rnaCalls))

    ##SNP pruning
    rnaCalls[rnaCalls == 0] <- NA
    
    maf <- apply(rnaCalls, 1, function(x) min(tabulate(x))/length(x))
    
    noninformative <- maf == 1 | !is.finite(maf) | maf < mafThr
    
    if(verbose) {
        hist(maf)
        message("noninformative RNA SNP calls: ", sum(noninformative))
    }
    rnaCalls <- rnaCalls[!noninformative,]
    
    rnaCalls  
}

DNACalls <- function(cohort, snps=NULL, DNAFile, type, verbose){

    suppressPackageStartupMessages({
        require(BBMRIomics)
        require(BiocParallel)
    })

    if(type == "HRC") {
        imputations <- getView("getImputations", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB, verbose=verbose)
        if(cohort != "ALL")
            ids <- unique(subset(imputations, imputation_reference == "HRC" & biobank_id == cohort)$imputation_id)
        else
            ids <- unique(subset(imputations, imputation_reference == "HRC")$imputation_id)
    }
    else if (type == "GoNL"){
        imputations <- getView("getIds", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB, verbose=verbose)
        if(cohort != "ALL")
            ids <- unique(subset(imputations, !is.na(gonl_id))$gonl_id)
        else
            ids <- unique(subset(imputations, biobank_id == cohort & !is.na(gonl_id))$gonl_id)
    }

    if(length(ids) == 0)
        stop("No imputed samples for this reference: ", type)

    if(is.null(snps)) {
        require(rtracklayer)
        snps <- import(DNAFile)
    }

    register(MulticoreParam(6))
    dnaCalls <- getGenotypes(ids, cohort, snps, type=type, geno="SM", BASE=VM_BASE_DATA)
}

relabelIntra <- function(x) {
    x <- gsub("original|has monozygotic twin|has repeated measurements|merged|replicate|rerun", "identical", x)
    x <- gsub("has child|has parent", "parentoffspring", x)
    x <- gsub("has dizygotic twin|has sib", "sibship", x)    
    as.character(x)
}

relabelInter <- function(x) {
    x <- gsub("original|has monozygotic twin|has repeated measurements|merged|replicate|rerun", "identical", x)
    x <- gsub("has dizygotic twin|has sib|has child|has parent", "unrelated", x)
    as.character(x)
}


getRelations <- function(type, verbose){

    ##obtain relations
    relations <- getView("getRelations", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB, verbose=verbose)
    relations <- subset(relations, !is.na(relation_type))

    if(length(unique(unlist(strsplit(type, "-")))) == 1) {
        if(type == "GoNL-GoNL") {

            relx <- getView("getImputations", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB, verbose=verbose)
            colnames(relx)[colnames(relx) == "gonl_id"] <- "idx"
            colnames(relations)[colnames(relations) == "gonl_id"] <- "idx"

        } else if (type == "HRC-HRC") {

            relx <- getView("getImputations", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB, verbose=verbose)
            relx <- subset(relx, imputation_reference == "HRC")
            colnames(relx)[colnames(relx) == "imputation_id"] <- "idx"

        } else if(type == "DNAm-DNAm") {

            relx <- getView("getMethylationRuns", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB, verbose=verbose)
            colnames(relx)[colnames(relx) == "run_id"] <- "idx"

        } else if (type == "RNA-RNA") {

            relx <- getView("getRNASeqRuns", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB, verbose=verbose)
            colnames(relx)[colnames(relx) == "run_id"] <- "idx"

        }

        relx <- relx[!is.na(relx$idx),]
        relx <- relx[!duplicated(relx$idx),]

        relations <- merge(relations, relx[,c("ids", "idx")], by="ids") #link idx
        relations <- relations[!is.na(relations$idx),]

        ##obtain reported relations
        bRels <- merge(relations, relx[, c("ids", "idx")], by.x="relation_id", by.y="ids")[,c("idx.x", "idx.y", "relation_type")]

        ##obtain technical relations
        tRels <- merge(relx, relx, by="ids")[, c("idx.x", "idx.y")]

    } else if (length(unique(unlist(strsplit(type, "-")))) == 2) {

        if( type == "DNAm-GoNL" | type == "GoNL-DNAm") {
            stop("not implemented yet!")

        } else if( type == "DNAm-HRC" | type == "HRC-DNAm") {

            relx <- getView("getMethylationRuns", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB, verbose=verbose)
            colnames(relx)[colnames(relx) == "run_id"] <- "idx"
            relx <- relx[!duplicated(relx$idx), ]

            rely <- getView("getImputations", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB, verbose=verbose)
            rely <- subset(rely, imputation_reference == "HRC")
            colnames(rely)[colnames(rely) == "imputation_id"] <- "idy"

        } else if( type == "RNA-HRC" | type == "HRC-RNA") {

            relx <- getView("getRNASeqRuns", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB, verbose=verbose)
            colnames(relx)[colnames(relx) == "run_id"] <- "idx"
            relx <- relx[!duplicated(relx$idx), ]

            rely <- getView("getImputations", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB, verbose=TRUE)
            rely <- subset(rely, imputation_reference == "HRC")
            colnames(rely)[colnames(rely) == "imputation_id"] <- "idy"

        } else if( type == "RNA-GoNL" | type == "GoNL-RNA") {
            stop("not implemented yet!")
        }

        relx <- relx[!is.na(relx$idx),]
        relx <- relx[!duplicated(relx$idx),]

        rely <- rely[!is.na(rely$idy),]
        rely <- rely[!duplicated(rely$idy),]

        bRelsx <- merge(relations, relx[,c("ids", "idx")], by="ids") #link idx
        bRelsx <- merge(bRelsx, rely[,c("ids", "idy")], by.x="relation_id", by.y="ids")[, c("idx", "idy", "relation_type")]

        bRelsy <- merge(relations, rely[,c("ids", "idy")], by="ids") #link idy
        bRelsy <- merge(bRelsy, relx[,c("ids", "idx")], by.x="relation_id", by.y="ids")[, c("idx", "idy", "relation_type")]

        bRels <- rbind(bRelsx, bRelsy)
        tRels <- merge(relx, rely, by="ids")[, c("idx", "idy")]
    }

    tRels$relation_type <- "identical"
    relations <- rbind(bRels, tRels)

    relations
}

genotyping <- function(typex, typey, filex, filey, cohort, out, verbose) {

    suppressPackageStartupMessages({
        require(BBMRIomics)
        require(GenomicRanges)
    })

    message("Start genotyping...")

    ##not all combinations are allowed (yet)
    type <- match.arg(paste(typex, typey, sep="-"),
                      choices = c("GoNL-GoNL", "HRC-HRC", "DNAm-DNAm", "RNA-RNA",
                                  "DNAm-HRC","DNAm-GoNL", "RNA-HRC","RNA-GoNL",
                                  "HRC-DNAm","GoNL-DNAm", "HRC-RNA","GoNL-RNA"))

    if(verbose) {
        message("Type comparison: ", type)
        message("Cohort: ", cohort)
        message("Verbose: ", verbose)
    }

    message("Extracting the data...")

    ##hard-coded
    ##RNAFile <- file.path(VM_BASE_ANALYSIS, "BBMRIomics/data/output.vcf")
    ##RNAFile <- file.path("~/output.vcf")
    ##DNAmFile <- file.path(VM_BASE_ANALYSIS, "BBMRIomics/data/DNAm_snps.RData")
    ##DNAFile <- file.path(VM_BASE_ANALYSIS, "BBMRIomics/data/final_list_50_SNPs.corrected.bed")

    if(any(typey %in% c("DNAm", "RNA")) & any(typex %in% c("HRC", "GoNL"))) {
        tmp <- typey
        typey <- typex
        typex <- tmp
    }

    ##read omic type(s)
    xCalls <- switch(typex,
                     "GoNL"= DNACalls(cohort, file = DNAFile, type="GoNL", verbose = verbose),
                     "HRC"= DNACalls(cohort, file = DNAFile, type="HRC", verbose = verbose),
                     "DNAm"= DNAmCalls(cohort, filex, verbose),
                     "RNA"= RNACalls(filex, verbose=verbose)) ##no subsetting on cohorts yet

    if(typey != typex) {
        if(typex == "RNA") {
            snps <- GRanges(seqnames=gsub(":.*$", "", rownames(xCalls)),
                            IRanges(start=as.integer(gsub("^.*:", "", rownames(xCalls))), width=1))
        }
        if(typex == "DNAm") {
            load(DNAmFile)
            snps <- DNAmSNPs[names(DNAmSNPs) %in% rownames(xCalls)]
            snpsC <- snpsG <- snps
            end(snpsC) <- start(snpsC)
            start(snpsG) <- end(snpsG)
            snps <- c(snpsC, snpsG)

            map <- names(snps)
            names(map) <- paste(seqnames(snps), start(snps), sep=":")
        }

        yCalls <- switch(typey,
                         "GoNL"= DNACalls(cohort, snps = snps, type="GoNL", verbose = verbose),
                         "HRC"= DNACalls(cohort, snps = snps, type="HRC", verbose = verbose),
                         "DNAm"= DNAmCalls(cohort, filey, verbose),
                         "RNA"= RNACalls(filey, verbose=verbose)) ##no subsetting on cohorts yet

        if(typex == "DNAm") {
            mid <- match(rownames(yCalls), names(map))
            rownames(yCalls) <- map[mid]
        }

    } else
        yCalls <- NULL

    message("Obtaining relations...")

    ##obtain relations
    type <- paste(typex, typey, sep="-")
    relations <- getRelations(type, verbose=verbose)

    if(typex == typey) {
        relations$relation_type <- relabelIntra(relations$relation_type)
        rHash <- hashRelations(relations, idx.col="idx.x", idy.col="idx.y")
    } else {
        relations$relation_type <- relabelInter(relations$relation_type)
        rHash <- hashRelations(relations, idx.col="idx", idy.col="idy")
    }

    message("Run allelesharing algorithm...")

    ##run allele sharing
    if(typex != typey & typex == "DNAm")
        data <- alleleSharing(x=xCalls, y=yCalls, rHash=rHash, phasing=TRUE, verbose=verbose)
    else
        data <- alleleSharing(x=xCalls, y=yCalls, rHash=rHash, verbose=verbose)
    
    if(verbose)
        print(head(data))

    message("Predict sample relations...")

    ##predict and output
    if(!is.null(out))
        fileName <- file.path(out, paste0("mismatches_", type, "_", cohort))

    if(!is.null(out)) {
        pdf(paste0(fileName, ".pdf"))
        mismatches <- predict(data)
        title(paste(cohort, "(", type, ")", collapse=""))
        dev.off()
    }
    else {
        mismatches <- predict(data)
        title(paste(cohort, "(", type, ")", collapse=""))
    }

    if(nrow(mismatches) > 0) {

        if(typex == "RNA")
            runs <- getView("getRNASeqRuns", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB, verbose=verbose)
        else if(typex == "DNAm")
            runs <- getView("getMethylationRuns", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB, verbose=verbose)

        mm <- merge(mismatches, runs[, c("run_id", "ids")], by.x="colnames.x", by.y="run_id")

        if(typey == typex)
            mm$ids.y <- runs$ids[match(mm$colnames.y, runs$run_id)]
        else if (typey == "HRC") {

            runs <- getView("getImputations", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB, verbose=verbose)
            runs <- subset(runs, imputation_reference == "HRC")
            colnames(runs)[colnames(runs) == "imputation_id"] <- "run_id"

        } else if (typey == "GoNL") {

            runs <- getView("getImputations", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB, verbose=verbose)
            colnames(runs)[colnames(runs) == "gonl_id"] <- "run_id"

        }

        mm <- merge(mm, runs[, c("run_id", "ids")], by.x="colnames.y", by.y="run_id")

        ##colnames(mm) <-  c("colnames.y", "colnames.x", "mean", "var", "relation", "predicted", "ids.x", "ids.y")
        mm <- mm[, c( "mean", "var", "relation", "predicted", "colnames.x", "ids.x", "colnames.y", "ids.y")]
        mm[, 1:2] <- round(mm[, 1:2], 3)
        if(!is.null(out))
            write.table(mm, file=paste0(fileName, ".txt"), row.names=FALSE, quote=FALSE, sep="\t")
        else
            return(invisible(mm))
    }

}
