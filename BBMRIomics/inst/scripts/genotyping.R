#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    require(optparse)
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

if (is.null(opt$filex) | !file.exists(opt$filex)){
    print_help(opt_parser)
    stop("At least one existing omic type file should be given.", call.=FALSE)
}

if (is.null(opt$typey) & !is.null(opt$filey)) {
    stop("You should provide to omic type of file:", opt$filey)
}

if (!is.null(opt$typey) & is.null(opt$filey)) {
    stop("You should provide a file for omic type:", opt$typey)
}

if (is.null(opt$typey)) {
    opt$typey <- opt$typex
} else {
    opt$typex <- match.arg(opt$typex, choices=c("GoNL", "HRC", "DNAm", "RNA"))
}


opt$cohort <- match.arg(opt$cohort, choices=c("ALL", "CODAM", "LL", "LLS", "NTR", "PAN", "RS"))

suppressPackageStartupMessages({
    require(BBBMRIomics)
    source(file.path(path.package("BBMRIomics"), "scripts/Genotyping_Helpers.R"), verbose=FALSE)
})

## typex <- "RNA"
## typey <- "RNA"
## filex <- "~/output.vcf"
## filey <- file.path(VM_BASE_ANALYSIS, "BBMRIomics/data/DNAm_snps.RData")
## cohort <- "ALL"
## verbose <- TRUE
## genotyping(typex, typey, filex, filey, cohort=cohort, out=NULL, verbose=verbose)

genotyping(typex=opt$typex, typey=opt$typey, filex=opt$filex, filey=opt$filey, cohort=opt$cohort, out=opt$out, verbose=opt$verbose)
