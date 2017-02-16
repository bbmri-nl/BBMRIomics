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
    data$relation <- 1L  ##unrelated == 1
    data$relation[mId] <- unlist(mget(pairs[mId], rHash, mode = "integer", ifnotfound=list(1L)), use.names=FALSE)  ##unrelated == 1    
    data$relation <- factor(data$relation, levels=c(1:4), labels=c("unrelated", "identical", "parentoffspring", "sibship"))
    invisible(droplevels(data))
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
    values <- as.integer(factor(relations[, rel.col])) ##unrelated == 1

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
