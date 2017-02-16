.onAttach <- function(libname, pkgname) {

    packageStartupMessage("##########################################################################################\n",
                          "##                                                                                      ##\n",                        
                          "## FACILITATING BBMRIOMICS DOWNSTREAM ANALYSES USING R                                  ##\n",
                          "##                                                                                      ##\n",
                          "##  DOCUMENTATION :                                                                     ##\n",
                          "##                  VIGNETTES : http://bios-vm.bbmrirp3-lumc.surf-hosted.nl/BBMRIomics/ ##\n",                          
                          "##                  WIKI      : http://wiki.bbmri.nl/wiki/BIOS_start-                   ##\n",
                          "##                                                                                      ##\n",
                          "##  QUESTIONS     :                                                                     ##\n",
                          "##                  FORUM : https://www.biostars.org/t/bbmri-nl-bios/                   ##\n",
                          "##                  EMAIL : m.van_iterson@lumc.nl or h.mei@lumc.nl                      ##\n",
                          "##                                                                                      ##\n",
                          "## ", pkgname, " (", packageVersion(pkgname), ")                                                                  ##\n",
                          "##########################################################################################\n")

    ##assign urls and directories
    assign("VM_BASE_DATA", "/virdir/Scratch/RP3_data", envir=as.environment(paste0("package:", pkgname)))
    assign("VM_BASE_ANALYSIS", "/virdir/Scratch/RP3_analysis", envir=as.environment(paste0("package:", pkgname)))    
    assign("SRM_BASE", "https://fly1.grid.sara.nl:2882/pnfs/grid.sara.nl/data/bbmri.nl/", envir=as.environment(paste0("package:", pkgname)))        
    assign("RP3_MDB", "https://metadatabase.bbmrirp3-lumc.surf-hosted.nl:6984/bios/", envir=as.environment(paste0("package:", pkgname)))
    assign("RP3_RDB", "https://metadatabase.bbmrirp3-lumc.surf-hosted.nl:6984/rp3_analysis/", envir=as.environment(paste0("package:", pkgname)))
    assign("RP4_DB", "https://db.metabolomicsdb-lumc.surf-hosted.nl/molgenis.R", envir=as.environment(paste0("package:", pkgname)))
    assign("RP3_BIOBANKS", c("RS", "PAN", "CODAM", "LLS", "LL", "NTR"), envir=as.environment(paste0("package:", pkgname)))

    ##TODO find a nicer way to show this information
    ##views/designs
    ## views <- list("freeze1RNASeq","freeze1Methylation", "freeze2RNASeq","freeze2Methylation",
    ##               "getFastq", "getIdat",
    ##               "getIds",
    ##               "allPhenotypes", "cellCounts", "minimalPhenotypes",
    ##               "getGenotypes", "getMethylationRuns", "getRNASeqRuns",
    ##               "rnaseqSamplesheet", "methylationSamplesheet",
    ##               "md5sum",
    ##               "getStats")

    ## assign("VIEWS", views, envir=globalenv())

    ##if exists assign usrpwd and proxy
    if(file.exists("~/.bbmriomics")) {
        ##drop these three lines in the next release
        usrpwdrp3 <- read.dcf("~/.bbmriomics", fields="usrpwd") ##for backward-compatibility                
        usrpwdrp3 <- read.dcf("~/.bbmriomics", fields="usrpwdrp3")
        usrpwdrp4 <- read.dcf("~/.bbmriomics", fields="usrpwdrp4")

        proxy <- read.dcf("~/.bbmriomics", fields="proxy")
        
        assign("RP3_MDB_USRPWD", as.character(usrpwdrp3), envir=as.environment(paste0("package:", pkgname)))
        assign("RP4_DB_USRPWD", as.character(usrpwdrp4), envir=as.environment(paste0("package:", pkgname)))
        assign("GRID_PROXY", as.character(proxy), envir=as.environment(paste0("package:", pkgname)))
    }
    else
        assign("RP3_MDB_USRPWD", "anonymous", envir=as.environment(paste0("package:", pkgname)))
}

##' Global variable
##'
##' @name RP3_MDB
##' @docType data
##' @author mvaniterson
##' @keywords data
RP3_MDB <- NULL

##' Global variable
##'
##' @name RP3_MDB_USRPWD
##' @docType data
##' @author mvaniterson
##' @keywords data
RP3_MDB_USRPWD <- NULL

##' Global variable
##'
##' @name RP4_DB
##' @docType data
##' @author mvaniterson
##' @keywords data
RP4_DB <- NULL

##' Global variable
##'
##' @name RP4_DB_USRPWD
##' @docType data
##' @author mvaniterson
##' @keywords data
RP4_DB_USRPWD <- NULL

##' Global variable
##'
##' @name VM_BASE_DATA
##' @docType data
##' @author mvaniterson
##' @keywords data
VM_BASE_DATA <- NULL

##' Global variable
##'
##' @name VM_BASE_ANALYSIS
##' @docType data
##' @author mvaniterson
##' @keywords data
VM_BASE_ANALYSIS <- NULL

##' Global variable
##'
##' @name SRM_BASE
##' @docType data
##' @author mvaniterson
##' @keywords data
SRM_BASE <- NULL

##' Global variable
##'
##' @name RP3_RDB
##' @docType data
##' @author mvaniterson
##' @keywords data
RP3_RDB <- NULL

##' Global variable
##'
##' @name GRID_PROXY
##' @docType data
##' @author mvaniterson
##' @keywords data
GRID_PROXY <- NULL

##' Global variable
##'
##' @name RP3_BIOBANKS
##' @docType data
##' @author mvaniterson
##' @keywords data
RP3_BIOBANKS <- NULL

molgenis.login <- function(...) NULL
molgenis.get <- function(...) NULL

