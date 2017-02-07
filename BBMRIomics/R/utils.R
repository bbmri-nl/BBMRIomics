.onLoad <- function(libname, pkgname) {

    packageStartupMessage("###################################################################################\n",
                          "## ", pkgname, " (", packageVersion(pkgname), ") package: Facilitating BIOS downstream analyses using R.    ##\n",
                          "##                                                                               ##\n",
                          "##  Documentation : http://wiki.bbmri.nl/wiki/BIOS_start-                        ##\n",
                          "##  Questions     : m.van_iterson@lumc.nl or h.mei@lumc.nl                       ##\n",
                          "###################################################################################\n")

    ##couchDB url
    couchHost <- "metadatabase.bbmrirp3-lumc.surf-hosted.nl"
    ##srm url
    srmUrl <- "fly1.grid.sara.nl:2882/pnfs/grid.sara.nl/data/"
    
    ##directories VM
    assign("RP3DATADIR", "/virdir/Scratch/RP3_data", envir=globalenv())

    ##summarize datasets
    assign("DATASETS",
           gsub("\\.RData", "", dir(file.path(find.package(pkgname), "data"), pattern="RData")),
           envir=globalenv())

    ##directories SRM
    assign("SRMBASE", paste0("https://", srmUrl, "bbmri.nl/"), envir=globalenv())    
    
    assign("MDB", paste0("https:", couchHost, ":6984/bios/"), envir=globalenv())
    assign("RDB", paste0("https:", couchHost, ":6984/rp3_analysis/"), envir=globalenv())

    ##views/designs
    views <- list("freeze1RNASeq","freeze1Methylation", "freeze2RNASeq","freeze2Methylation",
                  "getFastq", "getIdat",
                  "getIds",
                  "allPhenotypes", "cellCounts", "minimalPhenotypes",
                  "getGenotypes", "getMethylationRuns", "getRNASeqRuns",
                  "rnaseqSamplesheet", "methylationSamplesheet",
                  "md5sum",
                  "getStats")

    assign("VIEWS", views, envir=globalenv())

    ##if exists assign usrpwd/proxy
    if(file.exists("~/.biosrutils")) {
        usrpwdrp3 <- read.dcf("~/.biosrutils", fields="usrpwd") ##for backward-compatibility
        usrpwdrp3 <- read.dcf("~/.biosrutils", fields="usrpwdrp3")
        usrpwdrp4 <- read.dcf("~/.biosrutils", fields="usrpwdrp4")
        proxy <- read.dcf("~/.biosrutils", fields="proxy")
        assign("USRPWDRP3", usrpwdrp3, envir=globalenv())
        assign("USRPWDRP4", usrpwdrp4, envir=globalenv())
        assign("PROXY", proxy, envir=globalenv())
    }
    else
        assign("USRPWDRP3", "anonymous", envir=globalenv())
}

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

