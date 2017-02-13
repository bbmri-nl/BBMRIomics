
view <- function(viewname, db="http://127.0.0.1:5984/mdb_test/", converter="json2csv") {
    url <- file.path(db, "_design/couchdbapp/_list", converter, viewname)
    response <- RCurl::getURL(paste0(url, "?reduce=false"))
    tConn <- textConnection(response)
    on.exit(close(tConn))    
    invisible(read.table(tConn, sep=",", header=TRUE, na.string="null", as.is=TRUE))    
}

phenotypes <- view("phenotypes")

dim(phenotypes)
head(phenotypes)
str(phenotypes)


