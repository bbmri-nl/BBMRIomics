##' connect to a molgenis database using the R-API
##'
##' connect to a molgenis database using the R-API
##' @title connect to a molgenis database using the R-API
##' @param usrpwd username:password default reading from .biosrutils
##' @param url src source url molgenis database
##' @return setup molgenis connection
##' @author mvaniterson
##' @importFrom RCurl getURL
##' @export
molgenis.connect <- function(usrpwd, url){
    options(RCurlOptions = list(ssl.verifypeer=FALSE, ssl.verifyhost=FALSE))
    eval(expr = parse(text = getURL(url)), envir=globalenv())
    usrpwd <- unlist(strsplit(gsub("'", "", usrpwd), ":"))
    molgenis.login(usrpwd[1], usrpwd[2])
    message("\nRun 'ls()' to see the available functions to interact with the molgenis database!")
}

##' extract all data for a particulair entity from a molgenis database
##'
##' wrapper around molgenis.get
##' @title get all data
##' @param entity the name of the table to extract data from
##' @param nrows number of rows used in the chunking
##' @param verbose verbose option default TRUE
##' @return data.frame containing the whole table
##' @author mvaniterson
##' @export
molgenis.get.all <- function(entity, nrows=1000, verbose=TRUE){
    if(!exists("molgenis.get"))
        stop("First 'molgenis.connect' should be run successfully!")

    ##maybe fix this by adding ... to molgenis.get(as.is=TRUE)
    options(stringsAsFactors = FALSE)
    on.exit(options(stringsAsFactors=TRUE))

    data <- d <- molgenis.get(entity, start=0, num=nrows+1)
    start <- 1
    while(nrow(d) >= nrows) {
        start <- start + nrows
        d <- molgenis.get(entity, start=start, num=nrows)
        data <- rbind(data, d)
        if(verbose)
            message("Extracted ", nrow(data), " rows...")
    }
    data
}
