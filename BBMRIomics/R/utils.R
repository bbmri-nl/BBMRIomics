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
##' srm.path <- gsub("srm.*nl", SRM_BASE, srm.path) ##for curl access with need slightly different url##' 
##' srm.path
##' files
##' SRM2VM(srm.path, files, vm.path="~", proxy=GRID_PROXY)
##' ##get bam
##' bam <- getView("getBAM", usrpwd=RP3_MDB_USRPWD, url=RP3_RDB)
##' bam$path[1]
##' file <- basename(bam$path[1])
##' srm.path <- file.path(dirname(bam$path[1])) 
##' srm.path <- gsub("srm.*nl", SRM_BASE, srm.path) ##for curl access with need slightly different url##' 
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

