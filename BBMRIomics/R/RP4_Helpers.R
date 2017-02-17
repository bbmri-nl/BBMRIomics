##' connect to a molgenis database using the R-API
##'
##' connect to a molgenis database using the R-API
##' @title connect to a molgenis database using the R-API
##' @param usrpwd username:password default reading from .biosrutils
##' @param url src source url molgenis database
##' @param verbose show additional output
##' @return setup molgenis connection
##' @author mvaniterson
##' @importFrom httr POST config verbose content
##' @export
molgenis.connect <- function(usrpwd, url, verbose=FALSE){
    usrpwd <- unlist(strsplit(gsub("'", "", usrpwd), ":"))

    molgenis.api.url <<- url

    p <- if(verbose)
             POST(paste0(molgenis.api.url, "login/"),
                  body = list(username = usrpwd[1], password = usrpwd[2]),  encode = "json",
                  config(ssl_verifypeer = 0L, ssl_verifyhost = 0L),
                  verbose())
         else
             POST(paste0(molgenis.api.url, "login/"),
                  body = list(username = usrpwd[1], password = usrpwd[2]),  encode = "json",
                  config(ssl_verifypeer = 0L, ssl_verifyhost = 0L))

    molgenis.token <<-  content(p)$token
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
##' @importFrom RCurl getURL curlEscape
##' @export
molgenis.get.all <- function(entity, nrows=1000, verbose=TRUE){
   
    ##maybe fix this by adding ... to molgenis.get(as.is=TRUE)
    options(stringsAsFactors = FALSE)
    on.exit(options(stringsAsFactors=TRUE))

    data <- d <- .molgenis.get(entity, start=0, num=nrows+1)
    start <- 1
    while(nrow(d) >= nrows) {
        start <- start + nrows
        d <- .molgenis.get(entity, start=start, num=nrows)
        data <- rbind(data, d)
        if(verbose)
            message("Extracted ", nrow(data), " rows...")
    }
    
    ##convert "true"/"false" to TRUE/FALSE
    isBool <- function(x) all(unique(x) %in% c("false", "true", "", NA))
    
    cols <- which(apply(data, 2, isBool))    
    if(length(cols) > 1)
        data[, cols] <- apply(data[, cols], 2, as.logical)
    else if(length(cols) == 1)
        data[, cols] <- as.logical(data[, cols])
        
    data
}

.molgenis.get <- function(entity, q = NULL, start = 0, num = 1000, sortColumn= NULL, sortOrder = NULL, attributes = NULL) {
    options(RCurlOptions = list(ssl.verifypeer=FALSE, ssl.verifyhost=FALSE))
    
    url <- paste0(molgenis.api.url, "csv/", entity, "?molgenis-token=", molgenis.token, "&start=", start, "&num=", num, "&sortColumn=", sortColumn, "&sortOrder=", sortOrder)

    if (!is.null(q)) {
        url <- paste0(url, "&q=", curlEscape(q))
    }

    if (!is.null(attributes)) {
        url <- paste0(url, "&attributes=", curlEscape(paste0(attributes, collapse = ",")))
    }

    ## FIXME Check metadata for every column and set a colClass vector corresponding to the correct type
    ## EXAMPLE: column1 contains strings,
    ## characterClass <- c("character")
    ## names(characterClass) <- c("column1")
    ## read.csv(url, colClass = c(characterClass))

    csv <- getURL(url)
    dataFrame <- read.csv(textConnection(csv), comment.char="", na.strings=c(""," ","NA", -9))

    return (dataFrame)
}
