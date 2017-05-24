DNAmCalls <- function(cohort, verbose=FALSE, maxbatch=500){

    suppressPackageStartupMessages({
        require(BBMRIomics)
        require(minfi)
        require(DNAmArray)
        require(BiocParallel)
        require(FDb.InfiniumMethylation.hg19)
        require(omicsPrint)
    })

    samplesheets <- getView("methylationSamplesheet", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB, verbose=verbose)

    ##get location idat-files on VM
    path450k <- file.path(VM_BASE_DATA, "IlluminaHumanMethylation450k")
    samplesheets$biobank_id <- gsub("-.*$", "", samplesheets$ids)
    samplesheets$Basename <- with(samplesheets, file.path(path450k, "raw", Sentrix_Barcode,
                                                          paste(Sentrix_Barcode, Sentrix_Position, sep = "_")))

    samplesheets <- samplesheets[!duplicated(samplesheets$run_id),]

    runs <- getView("getMethylationRuns", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB, verbose=verbose)
    runs <- runs[runs$qc == "passed",]
    targets <- samplesheets[samplesheets$run_id %in% runs$run_id, ]

    if(cohort != "ALL")
        targets <- targets[targets$biobank_id == cohort,]

    data(hg19.GoNLsnps)
    cpgs <- unique(hg19.GoNLsnps$probe)

    if(nrow(targets) > maxbatch) {
        betas <- lapply(split(targets, 1+(1:nrow(targets))%/%maxbatch), function(targetsbatch) {
            RGset <- read.metharray.exp.par(targetsbatch, verbose=verbose)
            beta <- getBeta(RGset)
            rbind(beta[rownames(beta) %in% cpgs, ], getSnpBeta(RGset))
        })
        betas <- do.call('cbind', betas)
    } else {
        RGset <- read.metharray.exp.par(targets, verbose=verbose)
        betas <- getBeta(RGset)
        betas <- rbind(betas[rownames(betas) %in% cpgs,], getSnpBeta(RGset))
    }

    beta2genotype(betas)
}


RNACalls <- function(vcfFile, verbose){

    suppressPackageStartupMessages({
        require(VariantAnnotation)
        require(Rsamtools)
    })

    vcf <- readVcf(vcfFile, "hg19")
    rnaCalls <- genotypeToSnpMatrix(vcf)$genotypes
    rnaCalls <- t(matrix(as.numeric(rnaCalls), nrow=nrow(rnaCalls), ncol=ncol(rnaCalls), dimnames=dimnames(rnaCalls)))

    rownames(rnaCalls) <- gsub("_.*$", "", rownames(rnaCalls))
    colnames(rnaCalls) <- gsub("\\.variant.*$", "", colnames(rnaCalls))

    rnaCalls[rnaCalls == 0] <- NA

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
        if(cohort == "ALL")
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

    getGenotypes(ids, cohort, snps, type=type, geno="SM", BASE=VM_BASE_DATA)
}

relabelIntra <- function(x) {
    x <- gsub("original|has monozygotic twin|has repeated measurements|merged|replicate|rerun", "identical", x)
    ##x <- gsub("has child|has parent", "parentoffspring", x)
    ##x <- gsub("has dizygotic twin|has sib", "sibship", x)
    x <- gsub("inferred 1st degree family|2nd degree family|has child|has parent|has dizygotic twin|has sib", "family", x)    
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

            relx <- getView("getMethylationRuns", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB, verbose=verbose)
            colnames(relx)[colnames(relx) == "run_id"] <- "idx"
            relx <- relx[!duplicated(relx$idx), ]

            rely <- getView("getImputations", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB, verbose=TRUE)
            rely <- subset(rely, !is.na(gonl_id))
            colnames(rely)[colnames(rely) == "gonl_id"] <- "idy"

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

            relx <- getView("getRNASeqRuns", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB, verbose=verbose)
            colnames(relx)[colnames(relx) == "run_id"] <- "idx"
            relx <- relx[!duplicated(relx$idx), ]

            rely <- getView("getImputations", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB, verbose=TRUE)
            rely <- subset(rely, !is.na(gonl_id))
            colnames(rely)[colnames(rely) == "gonl_id"] <- "idy"

        } else if( type == "HRC-GoNL" | type == "GoNL-HRC") {

            relx <- getView("getImputations", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB, verbose=verbose)
            relx <- subset(relx, imputation_reference == "HRC")
            colnames(relx)[colnames(relx) == "imputation_id"] <- "idx"

            rely <- getView("getImputations", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB, verbose=TRUE)
            rely <- subset(rely, !is.na(gonl_id))
            colnames(rely)[colnames(rely) == "gonl_id"] <- "idy"

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
        require(omicsPrint)
    })

    message("Start genotyping...")

    ##not all combinations are allowed (yet)
    type <- match.arg(paste(typex, typey, sep="-"),
                      choices = c("GoNL-GoNL", "HRC-HRC", "DNAm-DNAm", "RNA-RNA",
                                  "DNAm-HRC","DNAm-GoNL", "RNA-HRC","RNA-GoNL",
                                  "HRC-DNAm","GoNL-DNAm", "HRC-RNA","GoNL-RNA",
                                  "GoNL-HRC", "HRC-GoNL"))

    if(verbose) {
        message("Type comparison: ", type)
        message("Cohort: ", cohort)
        message("Verbose: ", verbose)
    }

    message("Extracting the data...")

    if(any(typey %in% c("DNAm", "RNA")) & any(typex %in% c("HRC", "GoNL"))) {
        tmp <- typey
        typey <- typex
        typex <- tmp
    }

    if(typex %in% c("HRC", "GoNL"))
        snps <- rtracklayer::import(filex)

    ##read omic type(s)
    xCalls <- switch(typex,
                     "GoNL"= DNACalls(cohort, snps = snps, type="GoNL", verbose = verbose),
                     "HRC"= DNACalls(cohort, snps = snps, type="HRC", verbose = verbose),
                     "DNAm"= DNAmCalls(cohort, verbose),
                     "RNA"= RNACalls(filex, verbose=verbose)) ##no subsetting on cohorts yet

    if(typey != typex) {
        if(typex == "RNA") {
            snps <- GRanges(seqnames=gsub(":.*$", "", rownames(xCalls)),
                            IRanges(start=as.integer(gsub("^.*:", "", rownames(xCalls))), width=1))
        }
        if(typex == "DNAm") {
            data(hg19.GoNLsnps)
            snps <- hg19.GoNLsnps[hg19.GoNLsnps$probe %in% rownames(xCalls) &
                                  hg19.GoNLsnps$variantType == "SNP",]
            map <- snps$probe
            names(map) <- paste(snps$CHROM, snps$snpBeg, sep=":")
            names(map) <- gsub("chr", "", names(map))

            snps <- GRanges(seqnames = gsub(":.*$", "", names(map)),
                            IRanges(start=as.integer(gsub("^.*:", "", names(map))), width=1))
            snps <- sort(snps)
            snps <- snps[!duplicated(snps),]
        }

        yCalls <- switch(typey,
                         "GoNL"= DNACalls(cohort, snps = snps, type="GoNL", verbose = verbose),
                         "HRC"= DNACalls(cohort, snps = snps, type="HRC", verbose = verbose),
                         "DNAm"= DNAmCalls(cohort,  verbose),
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
        ##relations$relation_type <- relabelInter(relations$relation_type)
        relations$relation_type <- relabelIntra(relations$relation_type)
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
    if(!is.null(out)) {
        fileName <- file.path(out, paste0("mismatches_", type, "_", cohort))
    }

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

        if(!(typex %in% c("HRC", "GoNL"))) {

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

        }
        else {

            if (typex == "HRC") {

                runs <- getView("getImputations", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB, verbose=verbose)
                runs <- subset(runs, imputation_reference == "HRC")
                colnames(runs)[colnames(runs) == "imputation_id"] <- "run_id"

            } else if (typex == "GoNL") {

                runs <- getView("getImputations", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB, verbose=verbose)
                colnames(runs)[colnames(runs) == "gonl_id"] <- "run_id"

            }

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

        }

        mm <- mm[, c( "mean", "var", "relation", "predicted", "colnames.x", "ids.x", "colnames.y", "ids.y")]
        mm[, 1:2] <- round(mm[, 1:2], 3)
        mm <- mm[!duplicated(mm),]

        if(!is.null(out))
            write.table(mm, file=paste0(fileName, ".txt"), row.names=FALSE, quote=FALSE, sep="\t")
        else
            return(invisible(mm))
    }

}
