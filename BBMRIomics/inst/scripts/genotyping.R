#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    require(optparse)
})

option_list <- list(
    make_option(c("-x", "--typex"), type="character", default=NULL,
                help="Omic type of dataset 1 ('GoNL', 'HRC', 'DNAm', 'RNA')", metavar="character"),
    make_option(c("-y", "--typey"), type="character", default=NULL,
                help="Omic type of dataset 2 ('GoNL', 'HRC', 'DNAm', 'RNA')", metavar="character"),
    make_option(c("-c", "--cohort"), type="character", default="ALL",
                help="name cohort ('ALL', 'CODAM', 'LL', 'LLS', 'NTR', 'PAN', 'RS')", metavar="character"),
    make_option(c("-o", "--out"), type="character", default=".",
                help="output directory [default= %default]", metavar="character"),
    make_option(c("-v", "--verbose"), type="logical", default=FALSE,
                help="show additional output [default= %default]", metavar="charachter")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

##checking input
if (is.null(opt$typex)){
    print_help(opt_parser)
    stop("At least one omic type should be given.", call.=FALSE)
} else {
    opt$typex <- match.arg(opt$typex, choices=c("GoNL", "HRC", "DNAm", "RNA"))
}

if (is.null(opt$typey)) {
    opt$typey <- opt$typex
} else {
    opt$typex <- match.arg(opt$typex, choices=c("GoNL", "HRC", "DNAm", "RNA"))
}

opt$cohort <- match.arg(opt$cohort, choices=c("ALL", "CODAM", "LL", "LLS", "NTR", "PAN", "RS"))

suppressPackageStartupMessages({
    require(BIOSRutils)
    ##MDB <- "https://metadatabase.bbmrirp3-lumc.surf-hosted.nl:6984/bios"
    USRPWDRP3 <- "anonymous"
})

##some helper functions

DNAmCalls <- function(cohort, DNAmFile, verbose=FALSE, maxbatch=500){

    suppressPackageStartupMessages({
        require(BIOSRutils)
        require(minfi)
        require(Leiden450K)
        require(BiocParallel)
        require(FDb.InfiniumMethylation.hg19)
    })

    samplesheets <- getView("methylationSamplesheet", verbose=verbose)

    ##get location idat-files on VM
    path450k <- file.path(RP3DATADIR, "IlluminaHumanMethylation450k")
    samplesheets$biobank_id <- gsub("-.*$", "", samplesheets$ids)
    samplesheets$Basename <- with(samplesheets, file.path(path450k, "raw", Sentrix_Barcode,
                                                          paste(Sentrix_Barcode, Sentrix_Position, sep = "_")))

    samplesheets <- samplesheets[!duplicated(samplesheets$run_id),]

    runs <- getView("getMethylationRuns", verbose=verbose)
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
            RGset <- read.450k.exp.par(targetsbatch, verbose=verbose)
            beta <- getBeta(RGset)
            betas <- rbind(beta[rownames(beta) %in% cpgs, ], getSnpBeta(RGset))
        })
        betas <- do.call('cbind', betas)
    } else {
        RGset <- read.450k.exp.par(targets, verbose=verbose)
        betas <- getBeta(RGset)
        betas <- rbind(betas[rownames(betas) %in% cpgs,], getSnpBeta(RGset))
    }

    beta2genotype(betas)
}


RNACalls <- function(vcfFile){

    suppressPackageStartupMessages({
        require(VariantAnnotation)
        require(Rsamtools)
    })

    if(!file.exists(gsub("vcf.gz", "vcf.gz.tbi", vcfFile)))
        indexTabix(vcfFile, format="vcf")
    vcf <- readVcf(TabixFile(vcfFile), "hg19")
    rnaCalls <- genotypeToSnpMatrix(vcf)$genotypes
    rnaCalls <- t(matrix(as.numeric(rnaCalls), nrow=nrow(rnaCalls), ncol=ncol(rnaCalls), dimnames=dimnames(rnaCalls)))
    rownames(rnaCalls) <- gsub("_.*$", "", rownames(rnaCalls))
    rnaCalls
}

DNACalls <- function(cohort, snps=NULL, DNAFile, type, verbose){

    suppressPackageStartupMessages({
        require(BIOSRutils)
        require(BiocParallel)
    })

    if(type == "HRC") {
        imputations <- getView("getImputations", verbose=verbose)
        if(cohort != "ALL")
            ids <- unique(subset(imputations, imputation_reference == "HRC" & biobank_id == cohort)$imputation_id)
        else
            ids <- unique(subset(imputations, imputation_reference == "HRC")$imputation_id)
    }
    else if (type == "GoNL"){
        imputations <- getView("getIds", verbose=verbose)
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
    dnaCalls <- getGenotypes(ids, cohort, snps, type=type)
}


relabelIntra <- function(x) {
    x <- gsub("original|has monozygotic twin|has repeated measurements|merged|replicate|rerun", "identical", x)
    x <- gsub("has child|has parent", "parentoffspring", x)
    x <- gsub("has dizygotic twin|has sib", "sibship", x)
    factor(x, levels=c("unrelated", "identical", "parentoffspring", "sibship"))
}

relabelInter <- function(x) {
    x <- gsub("original|has monozygotic twin|has repeated measurements|merged|replicate|rerun", "identical", x)
    x <- gsub("has dizygotic twin|has sib|has child|has parent", "unrelated", x)
    factor(x, levels=c("unrelated", "identical"))
}


getRelations <- function(type, verbose){

    ##obtain relations
    relations <- getView("getRelations", verbose=verbose)
    relations <- subset(relations, !is.na(relation_type))

    if(length(unique(unlist(strsplit(type, "-")))) == 1) {
        if(type == "GoNL-GoNL") {

            relx <- getView("getImputations", verbose=verbose)
            colnames(relx)[colnames(relx) == "gonl_id"] <- "idx"
            colnames(relations)[colnames(relations) == "gonl_id"] <- "idx"

        } else if (type == "HRC-HRC") {

            relx <- getView("getImputations", verbose=verbose)
            relx <- subset(relx, imputation_reference == "HRC")
            colnames(relx)[colnames(relx) == "imputation_id"] <- "idx"

        } else if(type == "DNAm-DNAm") {

            relx <- getView("getMethylationRuns", verbose=verbose)
            colnames(relx)[colnames(relx) == "run_id"] <- "idx"

        } else if (type == "RNA-RNA") {

            relx <- getView("getRNASeqRuns", verbose=verbose)
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

            relx <- getView("getMethylationRuns", verbose=verbose)
            colnames(relx)[colnames(relx) == "run_id"] <- "idx"
            relx <- relx[!duplicated(relx$idx), ]

            rely <- getView("getImputations", verbose=verbose)
            rely <- subset(rely, imputation_reference == "HRC")
            colnames(rely)[colnames(rely) == "imputation_id"] <- "idy"

        } else if( type == "RNA-HRC" | type == "HRC-RNA") {

            relx <- getView("getRNASeqRuns", verbose=verbose)
            colnames(relx)[colnames(relx) == "run_id"] <- "idx"
            relx <- relx[!duplicated(relx$idx), ]

            rely <- getView("getImputations", verbose=TRUE)
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

genotyping <- function(typex, typey, cohort, out, verbose) {

    suppressPackageStartupMessages({
        require(BIOSRutils)
        require(GenomicRanges)
        source("/virdir/Backup/RP3_analysis/biosrutils/R/Genotyping.R", verbose=FALSE)
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
        message("Verbose option = ", verbose)
    }

    message("Extracting the data...")

    ##hard-coded
    RNAFile <- "/virdir/Scratch/RP3_analysis/rnasnp/combined.vcf.gz"
    DNAmFile <- "/virdir/Scratch/RP3_analysis/rnasnp/DNAm_snps.RData"
    DNAFile <- "/virdir/Scratch/RP3_analysis/rnasnp/SNP_positions.bed"

    if(any(typey %in% c("DNAm", "RNA")) & any(typex %in% c("HRC", "GoNL"))) {
        tmp <- typey
        typey <- typex
        typex <- tmp
    }

    ##read omic type(s)
    xCalls <- switch(typex,
                     "GoNL"= DNACalls(cohort, file = DNAFile, type="GoNL", verbose = verbose),
                     "HRC"= DNACalls(cohort, file = DNAFile, type="HRC", verbose = verbose),
                     "DNAm"= DNAmCalls(cohort, DNAmFile, verbose),
                     "RNA"= RNACalls(RNAFile)) ##no subsetting on cohorts yet

    if(typey != typex) {
        if(typex == "RNA") {
            snps <- GRanges(seqnames=gsub(":.*$", "", rownames(xCalls)),
                            IRanges(start=gsub("^.*:", "", rownames(xCalls)), width=1))
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
                         "DNAm"= DNAmCalls(cohort, DNAmFile, verbose),
                         "RNA"= RNACalls(RNAFile)) ##no subsetting on cohorts yet

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

    fileName <- file.path(out, paste0("mismatches_", type, "_", cohort))

    message("Predict sample relations...")

    ##predict and output
    pdf(paste0(fileName, ".pdf"))
    mismatches <- predict(data)
    title(paste(cohort, "(", type, ")", collapse=""))
    dev.off()

    if(nrow(mismatches) > 0) {

        if(typex == "RNA")
            runs <- getView("getRNASeqRuns", verbose=verbose)
        else if(typex == "DNAm")
            runs <- getView("getMethylationRuns", verbose=verbose)

        mm <- merge(mismatches, runs[, c("run_id", "ids")], by.x="colnames.x", by.y="run_id")

        if(typey == typex)
            mm$ids.y <- runs$ids[match(mm$colnames.y, runs$run_id)]
        else if (typey == "HRC") {
            
            runs <- getView("getImputations", verbose=verbose)
            runs <- subset(runs, imputation_reference == "HRC")
            colnames(runs)[colnames(runs) == "imputation_id"] <- "run_id"

        } else if (typey == "GoNL") {

            runs <- getView("getImputations", verbose=verbose)
            colnames(runs)[colnames(runs) == "gonl_id"] <- "run_id"

        }
        
        mm <- merge(mm, runs[, c("run_id", "ids")], by.x="colnames.y", by.y="run_id")
        
        ##colnames(mm) <-  c("colnames.y", "colnames.x", "mean", "var", "relation", "predicted", "ids.x", "ids.y")
        mm <- mm[, c( "mean", "var", "relation", "predicted", "colnames.x", "ids.x", "colnames.y", "ids.y")]
        mm[, 1:2] <- round(mm[, 1:2], 3)
        write.table(mm, file=paste0(fileName, ".txt"), row.names=FALSE, quote=FALSE, sep="\t")
    }

}

typex <- "DNAm"
typey <- "HRC"
cohort <- "NTR"
verbose <- TRUE

genotyping(typex=opt$typex, typey=opt$typey, cohort=opt$cohort, out=opt$out, verbose=opt$verbose)
