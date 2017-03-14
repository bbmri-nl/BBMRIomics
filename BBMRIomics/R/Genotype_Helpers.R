.scan2data.frame <- function(file, nrow, ncol, col.names=NULL, ...) {
    op <- options("stringsAsFactors" = FALSE)
    on.exit(op)
    df <- scan(file, nmax = nrow, ...)
    df <- lapply(df, type.convert, as.is=TRUE)
    attr(df, "row.names") <- .set_row_names(nrow)
    if(!is.null(col.names))
        names(df) <- col.names
    attr(df, "class") <- "data.frame"
    df
}

##' read dosages files impute2-transformed
##'
##' read dosages files impute2-transformed
##' @title read dosages files impute2-transformed
##' @param file character filename
##' @param yieldSize yieldSize for reading data in chunks
##' @param colClassesInfo describes the types of the columns
##' @param type data.frame GRanges or SummerizedExperiment
##' @param verbose default TRUE show progress message
##' @param ... additional arguments to scanTabix
##' @return data.frame
##' @author mvaniterson
##' @import SummarizedExperiment
##' @importFrom Rsamtools TabixFile headerTabix scanTabix yieldSize
##' @importFrom GenomicRanges GRanges
##' @importFrom IRanges IRanges
##' @importFrom S4Vectors DataFrame SimpleList
##' @importFrom utils type.convert
##' @export
##' @examples
##' \dontrun{
##' gzipped <- dir(file.path(VM_BASE_DATA, "GWAS_ImputationGoNLv5/dosages", RP3_BIOBANKS[1]),
##' pattern= "gz$", full.names=TRUE)
##' chunk <- read.dosages(gzipped[1], yieldSize=5000)
##' chunk[1:5, 1:10]
##' chunk <- read.dosages(gzipped[1], yieldSize=5000, type="GRanges")
##' chunk
##' chunk <- read.dosages(gzipped[1], yieldSize=5000, type="SummarizedExperiment")
##' chunk
##' colData(chunk)
##' rowRanges(chunk)
##' assay(chunk)[1:5, 1:5]
##' }
read.dosages <- function(file,  yieldSize=NULL, colClassesInfo = c("character", "character", "integer", "numeric",  "numeric",
                                                                   "numeric", "integer", "integer", "character", "integer", "character",
                                                                   "character"), type=c("data.frame", "GRanges", "SummarizedExperiment"), verbose=TRUE, ...)  {

    type <- match.arg(type)

    if(verbose)
        message("Reading chunk...")

    if(class(file) != "TabixFile")
        file <- TabixFile(file, yieldSize=yieldSize)

    header <- gsub("#", "", headerTabix(file)$header)
    header <- unlist(strsplit(header, "\\t"))

    value <- scanTabix(file, ...)[[1]]

    if(length(value) == 0)
        return(NULL)

    txt <- textConnection(value)
    on.exit(close(txt))

    if(type == "GRanges")
        colClasses <- c(colClassesInfo, rep("NULL", length(header)-length(colClassesInfo)))
    else
        colClasses <- c(colClassesInfo, rep("numeric", length(header)-length(colClassesInfo)))

    chunk <- .scan2data.frame(txt,
                              nrow=yieldSize(file),
                              ncol=length(header),
                              sep="\t",
                              what=as.list(colClasses),
                              quiet=TRUE,
                              skip=0,
                              col.names=header)

    switch(type,
           "data.frame"= chunk,
           "GRanges"= with(chunk,
                           GRanges(seqname = paste0("chr", chr),
                                   IRanges(as.integer(pos), width=1),
                                   rsid=rsid, ref=ref, alt=alt)),
           "SummarizedExperiment" = SummarizedExperiment(rowRanges = with(chunk,
                                                                          GRanges(seqname = paste0("chr", chr),
                                                                                  IRanges(as.integer(pos), width=1),
                                                                                  rsid=rsid, ref=ref, alt=alt)),
                                                         assays=SimpleList(dosage = data.matrix(chunk[, -c(1:length(colClassesInfo))])),
                                                         colData=DataFrame(gwas_id = colnames(chunk)[-c(1:length(colClassesInfo))]))
           )

}

.getHRC <- function(param, files, imputation_id, genotype){

    chr <- unique(as.character(seqnames(param)))
    fls <- grep(paste0("chr", chr, ".filtered.dose.vcf.gz"), files, value=TRUE)
    if(length(param) == 0 | is.null(fls))
        stop("No files found or no SNPs in input!")

    m <- lapply(fls, function(fl) { ##optionally multiple files per chromosome
        vcf <- readVcf(TabixFile(fl), "hg19", param=param)
        if(genotype=="SM") {
            m <- genotypeToSnpMatrix(vcf)$genotypes
            m <- t(matrix(as.numeric(m), nrow=nrow(m), ncol=ncol(m), dimnames=dimnames(m)))
        }
        else
            m <- geno(vcf)[[genotype]]
        m
    })
    m <- do.call("cbind", m)
    m <- m[, match(imputation_id, colnames(m)), drop=FALSE] ##return in proper order
    m
}

.getGONL <- function(param, files, imputation_id){

    chr <- unique(as.character(seqnames(param)))
    file <- grep(paste0("chr", chr, ".release5.raw_SNVs.vcf.gz"), files, value=TRUE)

    if(length(param) == 0 | is.null(file)) return(NULL)

    vcf <- readVcf(TabixFile(file), "hg19", param=param)
    m <- genotypeToSnpMatrix(vcf)$genotypes
    m <- t(matrix(as.numeric(m), nrow=nrow(m), ncol=ncol(m), dimnames=dimnames(m)))
    rownames(m) <- paste(seqnames(vcf@rowRanges), start(vcf@rowRanges), sep=":")
    m[, match(imputation_id, colnames(m)), drop=FALSE] ##return in proper order
}

##' extract genotypes from vcf-files
##'
##' extract genotypes from vcf-files
##' given selection of SNPs and samples
##' @title extract genotypes from vcf-files
##' @param imputation_id imputation identifier
##' @param biobank biobank_id
##' @param snps GRanges with snps
##' @param type imputation type either "GoNL", "HRC" or "GoNLv5"
##' @param geno extract either genotypes, dosages, genotype likelihoods or as snpMatrix
##' @param BASE genotype data location default e.g. VM_BASE_DATA
##' @param ... optional BPPARAM arguments
##' @return matrix with genotypes
##' @author mvaniterson
##' @importFrom VariantAnnotation readVcf genotypeToSnpMatrix geno
##' @importFrom Rsamtools TabixFile
##' @importFrom BiocParallel bplapply
##' @importFrom GenomicRanges split
##' @importFrom GenomeInfoDb mapSeqlevels
##' @export
getGenotypes <- function(imputation_id, biobank=c("ALL", "CODAM", "LL", "LLS", "NTR", "RS", "PAN"), snps, type=c("GoNL", "HRC", "GoNLv5"), geno=c("GT", "DS", "GP", "SM"), BASE, ...){

    type <- match.arg(type)
    biobank <- match.arg(biobank)
    geno <- match.arg(geno)

    ##snps should be of type GRanges
    seqlevels(snps) <- mapSeqlevels(seqlevels(snps), "NCBI")
    snps <- split(snps, as.character(seqnames(snps)))

    if(type == "HRC") {

        if(biobank == "ALL")
            vcfs <- dir(file.path(BASE, "HRC_Imputation"), pattern="filtered.dose.vcf.gz$", full.names=TRUE, recursive=TRUE)
        else
            vcfs <- dir(file.path(BASE, "HRC_Imputation", biobank), pattern="filtered.dose.vcf.gz$", full.names=TRUE, recursive=TRUE)


        ##for(fl in vcfs) indexTabix(fl, format="vcf") ##if vcf are not indexed!
        ##TODO Bioconductor devel (bioc-3.4/R-3.3.0) contains `GenomicFiles` with vcfstack a nicer solution?

        if(length(snps) > 1) {
            genotypes <- bplapply(snps, .getHRC, files=vcfs, imputation_id = as.character(imputation_id), genotype = geno, ...)
            genotypes <- do.call("rbind", genotypes)
        } else {
            genotypes <- .getHRC(snps[[1]], files=vcfs, imputation_id = as.character(imputation_id), genotype = geno, ...)
        }

    } else if(type == "GoNL") {
        vcfs <- dir(file.path(BASE, "gonl-snv-release-5.4"), pattern=".vcf.gz$", full.names=TRUE, recursive=TRUE)
        genotypes <- bplapply(snps, .getGONL, files=vcfs, imputation_id = as.character(imputation_id), genotype = geno)
    }
    else if(type == "GoNLv5")
        stop("Not implemented yet!")

    genotypes
}
