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
                          "##                  EMAIL : h.mei@lumc.nl, d.cats@lumc.nl                               ##\n",
                          "##                                                                                      ##\n",
                          "## ", pkgname, " (", packageVersion(pkgname), ")                                                                   ##\n",
                          "##########################################################################################\n")

    ##assign urls and directories
    configFile <- file.path(libname, "BBMRIomics/configure", "bbmriomics.conf")

    if(!file.exists(configFile)) {
        configFile <- file.path(libname, "BBMRIomics/inst/configure", "bbmriomics.conf")
        if(!file.exists(configFile))
            stop(configFile, " doesn't exists!")
    }

    cfs <- read.dcf(configFile,
                    fields=c("VM_BASE_DATA", "VM_BASE_ANALYSIS", "SRM_BASE", "RP3_MDB", "RP3_RDB", "RP4_DB", "SQL_DB"))[1,]

    tmp <- mapply(FUN = function(key, value)
        assign(key, as.character(value), envir = as.environment(paste0("package:", pkgname))),
        key = names(cfs), value = cfs)

    ##assuming this wil not change
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
    configFiles <- c("~/.bbmriomics", "~/.biosrutils") ##for backward compatibility BIOSRutils
    if(any(file.exists(configFiles))) {
        configFile <- configFiles[file.exists(configFiles)]
        usrpwdrp3 <- read.dcf(configFile, fields="usrpwd") ##for backward compatibility
        usrpwdrp3 <- read.dcf(configFile, fields="usrpwdrp3")
        usrpwdrp4 <- read.dcf(configFile, fields="usrpwdrp4")

        proxy <- read.dcf(configFile, fields="proxy")

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
##' @name SQL_DB
##' @docType data
##' @author dcats
##' @keywords data
SQL_DB <- NULL

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
