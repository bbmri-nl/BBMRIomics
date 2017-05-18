#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    require(optparse)
    require(BiocParallel)
    require(BBMRIomics)
    ##source(file.path(path.package("BBMRIomics"), "scripts/Genotyping_Helpers.R"), verbose=FALSE)
    source(file.path("/virdir/Backup/RP3_analysis/BBMRIomics/BBMRIomics/inst", "scripts/Genotyping_Helpers.R"), verbose=FALSE)
})

option_list <- list(
    make_option(c("-x", "--typex"), type="character", default=NULL,
                help="Omic type of dataset 1 ('GoNL', 'HRC', 'DNAm', 'RNA')", metavar="character"),
    make_option(c("-y", "--typey"), type="character", default=NULL,
                help="Omic type of dataset 2 ('GoNL', 'HRC', 'DNAm', 'RNA')", metavar="character"),
    make_option(c("-fx", "--filex"), type="character", default=NULL,
                help="File (full path) for typex", metavar="character"),
    make_option(c("-fy", "--filey"), type="character", default=NULL,
                help="File (full path) for typey", metavar="character"),
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

if (!is.null(opt$filex)) {
    if(!file.exists(opt$filex)){
        print_help(opt_parser)
        stop("filex does not exists!", call.=FALSE)
    }
}

if (!is.null(opt$filey)) {
    if(!file.exists(opt$filey)){
        print_help(opt_parser)
        stop("filey does not exists!", call.=FALSE)
    }
}

if (is.null(opt$typey)) {
    opt$typey <- opt$typex
} else {
    opt$typey <- match.arg(opt$typey, choices=c("GoNL", "HRC", "DNAm", "RNA"))
}

opt$cohort <- match.arg(opt$cohort, choices=c("ALL", "CODAM", "LL", "LLS", "NTR", "PAN", "RS"))

suppressPackageStartupMessages({  
    source(file.path("/virdir/Backup/RP3_analysis/BBMRIomics/BBMRIomics/inst", "scripts/Genotyping_Helpers.R"), verbose=FALSE)  
})

## opt <- list()
## opt$typex <- "GoNL"
## opt$typey <- "HRC"
## opt$filex <- "/virdir/Scratch/RP3_analysis/SwapDetection/HighQualPositions.GCRh37.bed"
## opt$filey <- NULL
## opt$cohort <- "ALL"
## opt$verbose <- TRUE
## opt$out <- "/virdir/Backup/RP3_analysis/SwapDetection/"

library(BiocParallel)
register(MulticoreParam(22))

genotyping(typex=opt$typex, typey=opt$typey, filex=opt$filex, filey=opt$filey, cohort=opt$cohort, out=opt$out, verbose=opt$verbose)

## typex <- "HRC"
## typey <- "GoNL"
## filex <- "/virdir/Scratch/RP3_analysis/SwapDetection/HighQualPositions.GCRh37.bed"
## filey <- NULL
## cohort <- "ALL"
## verbose <- TRUE
## out <- "/virdir/Backup/RP3_analysis/SwapDetection/"
