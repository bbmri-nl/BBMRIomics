##' Utility function to copy files from SRM to VM
##'
##' details follow
##' @title copy file from SRM2VM
##' @param srm.path srm path
##' @param files file or vector of files
##' @param vm.path vm path
##' @param proxy a valid proxy from example PROXY from your global
##' envirnoment
##' @return curl stdout message
##' @author mvaniterson
##' @export
##' @examples
##' \dontrun{
##' ##get fastq
##' fastq <- getView("getFastq", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB)
##' file <- basename(fastq$R1[1])
##' srm.path <- file.path(dirname(fastq$R1[1])) 
##' srm.path <- gsub("srm.*nl", SRM_BASE, srm.path) ##for curl access with need slightly different url 
##' srm.path
##' files
##' SRM2VM(srm.path, files, vm.path="~", proxy=GRID_PROXY)
##' ##get bam
##' bam <- getView("getBAM", usrpwd=RP3_MDB_USRPWD, url=RP3_RDB)
##' bam$path[1]
##' file <- basename(bam$path[1])
##' srm.path <- dirname(bam$path[1])
##' srm.path <- gsub("srm.*nl", SRM_BASE, srm.path) ##for curl access with need slightly different url 
##' srm.path
##' file
##' SRM2VM(srm.path, file, vm.path="~", proxy=GRID_PROXY)
##' }
SRM2VM <- function(srm.path, files, vm.path, proxy) {
    message(paste("\n############################################################################\n",
                  "## Besure the grid proxy is up-to-date!                                     ##\n",
                  "## Otherwise generate grid-proxy: `startGridSession bbmri.nl:/bbmri.nl/RP3` ##\n",
                  "## and copy proxy to VM by type this command from VM:                       ##\n",
                  "## scp `username@ui.grid.sara.nl:/tmp/yourproxy /tmp`                       ##\n",
                  "##############################################################################\n"))

    ##check file names using srmls  srm://srm.grid.sara.nl/pnfs ...
    for(file in files)
    {
        if(!file.exists(file.path(vm.path, dirname(file))))
            dir.create(file.path(vm.path, dirname(file)))

        cmd <- paste0("curl -k --CApath /etc/grid-security/certificates/ -E ", proxy, " -L ", file.path(srm.path, file), " > ", file.path(vm.path, file))
        print(cmd)
        system(cmd)
    }

    message(paste("Files copy to:", paste(file.path(vm.path, files)), collapse=""))
    invisible(paste(file.path(vm.path, files)))
}

##' head for matrices and data.frames as well
##'
##' optional p arguments selectes the number of columns to show
##' @title head
##' @param ... usual arguments to head
##' @param p numebr of columns to show
##' @return selection of rows and columns (in case matrix or data.frame)
##' @author mvaniterson
##' @importFrom utils head
##' @export
head <- function(..., p=6) {
    h <- utils::head(...)
    if(!is.null(dim(h))){
        p <- min(ncol(h), p)
        if(p < ncol(h))
            message("'head' showing ", p, " out of ", ncol(h), " columns!")
        h <- h[,1:p]
    }
    h
}


#' Load BBMRI data from the research drive
#'
#' @param dataset The dataset to be loaded.
#' @param envir The environment the dataset will be loaded into.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # list the available datasets
#' data(package="BBMRIomics")
#'
#' # load in a dataset
#' bbmri.data(metabolomics_RP3RP4_overlap)
#' }
bbmri.data <- function(dataset, envir = parent.frame(2)) {
    if (!(is.null(dataset) || typeof(dataset) == "character")) {
        stop("dataset should be null or character")
    }
    if (is.null(dataset)) {
        data.name = deparse(substitute(metabolomics_RP3RP4_overlap))
    } else {
        data.name = dataset
    }
    paths <- list(
        metabolomics_RP3RP4_overlap="RP4/metabolomics_RP3RP4_overlap.RData",
        methData_Betas_BIOS_F2_cleaned="IlluminaHumanMethylation450k/450k/methData_Betas_BIOS_F2_cleaned.RData",
        methData_Betas_BIOS_Freeze2_unrelated="IlluminaHumanMethylation450k/450k/methData_Betas_BIOS_Freeze2_unrelated.RData",
        methData_Betas_CODAM_F2_cleaned="IlluminaHumanMethylation450k/450k/CODAM/methData_Betas_CODAM_F2_cleaned.RData",
        methData_Betas_CODAM_Freeze2_unrelated="IlluminaHumanMethylation450k/450k/CODAM/methData_Betas_CODAM_Freeze2_unrelated.RData",
        methData_Betas_LLS_F2_cleaned="IlluminaHumanMethylation450k/450k/LLS/methData_Betas_LLS_F2_cleaned.RData",
        methData_Betas_LLS_Freeze2_unrelated="IlluminaHumanMethylation450k/450k/LLS/methData_Betas_LLS_Freeze2_unrelated.RData",
        methData_Betas_LL_F2_cleaned="IlluminaHumanMethylation450k/450k/LL/methData_Betas_LL_F2_cleaned.RData",
        methData_Betas_LL_Freeze2_unrelated="IlluminaHumanMethylation450k/450k/LL/methData_Betas_LL_Freeze2_unrelated.RData",
        methData_Betas_NTR_F2_cleaned="IlluminaHumanMethylation450k/450k/NTR/methData_Betas_NTR_F2_cleaned.RData",
        methData_Betas_NTR_Freeze2_unrelated="IlluminaHumanMethylation450k/450k/NTR/methData_Betas_NTR_Freeze2_unrelated.RData",
        methData_Betas_PAN_F2_cleaned="IlluminaHumanMethylation450k/450k/PAN/methData_Betas_PAN_F2_cleaned.RData",
        methData_Betas_PAN_Freeze2_unrelated="IlluminaHumanMethylation450k/450k/PAN/methData_Betas_PAN_Freeze2_unrelated.RData",
        methData_Betas_RS_F2_cleaned="IlluminaHumanMethylation450k/450k/RS/methData_Betas_RS_F2_cleaned.RData",
        methData_Betas_RS_Freeze2_unrelated="IlluminaHumanMethylation450k/450k/RS/methData_Betas_RS_Freeze2_unrelated.RData",
        methData_Mvalues_BIOS_F2_cleaned="IlluminaHumanMethylation450k/450k/methData_Mvalues_BIOS_F2_cleaned.RData",
        methData_Mvalues_BIOS_Freeze2_unrelated="IlluminaHumanMethylation450k/450k/methData_Mvalues_BIOS_Freeze2_unrelated.RData",
        methData_Mvalues_CODAM_F2_cleaned="IlluminaHumanMethylation450k/450k/CODAM/methData_Mvalues_CODAM_F2_cleaned.RData",
        methData_Mvalues_CODAM_Freeze2_unrelated="IlluminaHumanMethylation450k/450k/CODAM/methData_Mvalues_CODAM_Freeze2_unrelated.RData",
        methData_Mvalues_LLS_F2_cleaned="IlluminaHumanMethylation450k/450k/LLS/methData_Mvalues_LLS_F2_cleaned.RData",
        methData_Mvalues_LLS_Freeze2_unrelated="IlluminaHumanMethylation450k/450k/LLS/methData_Mvalues_LLS_Freeze2_unrelated.RData",
        methData_Mvalues_LL_F2_cleaned="IlluminaHumanMethylation450k/450k/LL/methData_Mvalues_LL_F2_cleaned.RData",
        methData_Mvalues_LL_Freeze2_unrelated="IlluminaHumanMethylation450k/450k/LL/methData_Mvalues_LL_Freeze2_unrelated.RData",
        methData_Mvalues_NTR_F2_cleaned="IlluminaHumanMethylation450k/450k/NTR/methData_Mvalues_NTR_F2_cleaned.RData",
        methData_Mvalues_NTR_Freeze2_unrelated="IlluminaHumanMethylation450k/450k/NTR/methData_Mvalues_NTR_Freeze2_unrelated.RData",
        methData_Mvalues_PAN_F2_cleaned="IlluminaHumanMethylation450k/450k/PAN/methData_Mvalues_PAN_F2_cleaned.RData",
        methData_Mvalues_PAN_Freeze2_unrelated="IlluminaHumanMethylation450k/450k/PAN/methData_Mvalues_PAN_Freeze2_unrelated.RData",
        methData_Mvalues_RS_F2_cleaned="IlluminaHumanMethylation450k/450k/RS/methData_Mvalues_RS_F2_cleaned.RData",
        methData_Mvalues_RS_Freeze2_unrelated="IlluminaHumanMethylation450k/450k/RS/methData_Mvalues_RS_Freeze2_unrelated.RData",
        rnaSeqData_ReadCounts_BIOS_Freeze2_GoNL="RNASeq/v2.1.3/gene_read/rnaSeqData_ReadCounts_BIOS_Freeze2_GoNL.RData",
        rnaSeqData_ReadCounts_BIOS_Freeze2_GoNL_GRCh38="RNASeq/GRCh38/gene_read/rnaSeqData_ReadCounts_BIOS_Freeze2_GoNL_GRCh38.RData",
        rnaSeqData_ReadCounts_BIOS_Freeze2_unrelated="RNASeq/v2.1.3/gene_read/rnaSeqData_ReadCounts_BIOS_Freeze2_unrelated.RData",
        rnaSeqData_ReadCounts_BIOS_Freeze2_unrelated_GRCh38="RNASeq/GRCh38/gene_read/rnaSeqData_ReadCounts_BIOS_Freeze2_unrelated_GRCh38.RData",
        rnaSeqData_ReadCounts_BIOS_cleaned="RNASeq/v2.1.3/gene_read/rnaSeqData_ReadCounts_BIOS_cleaned.RData",
        rnaSeqData_ReadCounts_CODAM_Freeze2_unrelated="RNASeq/v2.1.3/gene_read/rnaSeqData_ReadCounts_CODAM_Freeze2_unrelated.RData",
        rnaSeqData_ReadCounts_CODAM_Freeze2_unrelated_GRCh38="RNASeq/GRCh38/gene_read/rnaSeqData_ReadCounts_CODAM_Freeze2_unrelated_GRCh38.RData",
        rnaSeqData_ReadCounts_CODAM_cleaned="RNASeq/v2.1.3/gene_read/rnaSeqData_ReadCounts_CODAM_cleaned.RData",
        rnaSeqData_ReadCounts_GoNL="RNASeq/v2.1.3/gene_read/rnaSeqData_ReadCounts_GoNL.RData",
        rnaSeqData_ReadCounts_LLS_Freeze2_unrelated="RNASeq/v2.1.3/gene_read/rnaSeqData_ReadCounts_LLS_Freeze2_unrelated.RData",
        rnaSeqData_ReadCounts_LLS_Freeze2_unrelated_GRCh38="RNASeq/GRCh38/gene_read/rnaSeqData_ReadCounts_LLS_Freeze2_unrelated_GRCh38.RData",
        rnaSeqData_ReadCounts_LLS_cleaned="RNASeq/v2.1.3/gene_read/rnaSeqData_ReadCounts_LLS_cleaned.RData",
        rnaSeqData_ReadCounts_LL_Freeze2_unrelated="RNASeq/v2.1.3/gene_read/rnaSeqData_ReadCounts_LL_Freeze2_unrelated.RData",
        rnaSeqData_ReadCounts_LL_Freeze2_unrelated_GRCh38="RNASeq/GRCh38/gene_read/rnaSeqData_ReadCounts_LL_Freeze2_unrelated_GRCh38.RData",
        rnaSeqData_ReadCounts_LL_cleaned="RNASeq/v2.1.3/gene_read/rnaSeqData_ReadCounts_LL_cleaned.RData",
        rnaSeqData_ReadCounts_NTR_Freeze2_unrelated="RNASeq/v2.1.3/gene_read/rnaSeqData_ReadCounts_NTR_Freeze2_unrelated.RData",
        rnaSeqData_ReadCounts_NTR_Freeze2_unrelated_GRCh38="RNASeq/GRCh38/gene_read/rnaSeqData_ReadCounts_NTR_Freeze2_unrelated_GRCh38.RData",
        rnaSeqData_ReadCounts_NTR_cleaned="RNASeq/v2.1.3/gene_read/rnaSeqData_ReadCounts_NTR_cleaned.RData",
        rnaSeqData_ReadCounts_PAN_Freeze2_unrelated="RNASeq/v2.1.3/gene_read/rnaSeqData_ReadCounts_PAN_Freeze2_unrelated.RData",
        rnaSeqData_ReadCounts_PAN_Freeze2_unrelated_GRCh38="RNASeq/GRCh38/gene_read/rnaSeqData_ReadCounts_PAN_Freeze2_unrelated_GRCh38.RData",
        rnaSeqData_ReadCounts_PAN_cleaned="RNASeq/v2.1.3/gene_read/rnaSeqData_ReadCounts_PAN_cleaned.RData",
        rnaSeqData_ReadCounts_RS_Freeze2_unrelated="RNASeq/v2.1.3/gene_read/rnaSeqData_ReadCounts_RS_Freeze2_unrelated.RData",
        rnaSeqData_ReadCounts_RS_Freeze2_unrelated_GRCh38="RNASeq/GRCh38/gene_read/rnaSeqData_ReadCounts_RS_Freeze2_unrelated_GRCh38.RData",
        rnaSeqData_ReadCounts_RS_cleaned="RNASeq/v2.1.3/gene_read/rnaSeqData_ReadCounts_RS_cleaned.RData"
    )
    path <- file.path("~/researchdrive/RSC BIOS/RP3_data", paths[[data.name]])
    if (length(path) == 0) {
        stop("unknown dataset")
    }
    if (file.exists(path)) {
        load(path, envir = envir)
    } else {
        stop("dataset file does not exist")
    }
}

