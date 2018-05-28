##' get view from metadatabase
##'
##' @title get view from metadatabase
##' @param view Name of the to be retrieved view.
##' @param url Unique resource location of the metadatabase.
##' @param usrpwd Username and password concatenated by colon 'usr:pwd'.
##' This defaults to 'anonymous', which will cause a locally stored versions 
##' of the view to be retrieved.
##' @param port The port to be used to connect to the database. Defaults to 
##' 5432.
##' @param db The database to be connected to. Defaults to rp3_rp4_meta.
##' @importFrom RPostgreSQL dbDriver dbConnect dbGetQuery dbDisconnect
##' @export
##' @return data.frame containing the view
##' @author dcats
##' @examples
##' \dontrun{ 
##' #This example will not work yet, as RP3_MDB does not point at the SQL MDB.
##' view <- getSQLview("minimalphenotypes", url=RP3_MDB, usrpwd=RP3_MDB_USRPWD)
##' }
getSQLview <- function(view, url, usrpwd="anonymous", port=5432, 
                       db="rp3_rp4_meta") {
    views <- c("getfastq", "persontogwas_includingmztwins", "getidat",
               "freeze1rnaseq", "freeze2rnaseq", "freeze1methylation",
               "freeze2methylation", "getids", "allphenotypes", "cellcounts",
               "minimalphenotypes", "getimputatins", "getmethylationruns",
               "getrnaseqruns", "methylationsamplesheet", "rnaseqsamplesheet",
               "getrelations")
    
    if (! tolower(view) %in% views){
        stop(view, " not known!")
    }
    
    if(usrpwd != "anonymous") {
        usrpwd <- strsplit(usrpwd, ":")[[1]]
        drv <- dbDriver("PostgreSQL")
        con <- dbConnect(drv, dbname=db, host=url, port=port, user=usrpwd[1], 
                         password=usrpwd[2])
        tbl <- dbGetQuery(con, paste0("SELECT * FROM ", view))
        dbDisconnect(con)
        return(tbl)
    } else {
        message("No username and password provided for the MDB, using stored views!")
        
        pckg.path <- path.package("BBMRIomics")
        #pckg.path <- "~/testDataFolder"
        
        load(paste0(pckg.path, "/", tolower(view), ".Rdata"))
        return(tbl)
    }
}


##' update locally stored views
##'
##' Retrieve the views from the metadatabase and store them locally
##' as data.frames in Rdata files. Files will be stored in the 
##' BBMRIomics installation directory.
##' @title update locally stored views
##' @param url The url on which the metadatabase is hosted.
##' @param usrpwd The metadatabase username and password. Formatted as 
##' <username>:<password>.
##' @param port The port though which the database can be accessed. 
##' Defaults to 5432.
##' @param db Name of the database from which the views will be retrieved.
##' Defaults to "rp3_rp4_meta".
##' @author dcats
##' @importFrom RPostgreSQL dbDriver dbConnect dbGetQuery dbDisconnect
updateViews <- function(url, usrpwd, port=5432, db="rp3_rp4_meta") {
    views <- c("getfastq", "persontogwas_includingmztwins", "getidat",
               "freeze1rnaseq", "freeze2rnaseq", "freeze1methylation",
               "freeze2methylation", "getids", "allphenotypes", "cellcounts",
               "minimalphenotypes", "getimputatins", "getmethylationruns",
               "getrnaseqruns", "methylationsamplesheet", "rnaseqsamplesheet",
               "getrelations") #TODO make this automatic (and store list somewhere)?
    
    usrpwd <- strsplit(usrpwd, ":")[[1]]
    
    drv <- dbDriver("PostgreSQL")
    con <- dbConnect(drv, dbname=db, host=url, port=port, user=usrpwd[1],
                     password=usrpwd[2])
    
    pckg.path <- path.package("BBMRIomics")
    #pckg.path <- "~/testDataFolder"
    
    for (view in views){
        tbl <- dbGetQuery(con, paste0("SELECT * FROM ", view))
        save(tbl, file=paste0(pckg.path, "/", view, ".Rdata"))
    }
    
    dbDisconnect(con)
    return()
}

##' send a query to the database
##'
##' @title send a query to the database
##' @param query The query to be send to the database.
##' @param url The url on which the database is hosted.
##' @param usrpwd The database username and password. Formatted as 
##' <username>:<password>.
##' @param port The port though which the database can be accessed. 
##' Defaults to 5432.
##' @param db Name of the database from which the views will be retrieved.
##' Defaults to "rp3_rp4_meta".
##' @author dcats
##' @importFrom RPostgreSQL dbDriver dbConnect dbGetQuery dbDisconnect
##' @examples 
##' \dontrun{
##' visits <- .runQuery("SELECT * FROM visit;", RP3_MDB, RP3_MDB_USRPWD)
##' }
.runQuery <- function(query, url, usrpwd, port=5432, db="rp3_rp4_meta"){
    usrpwd <- strsplit(usrpwd, ":")[[1]]
    
    drv <- dbDriver("PostgreSQL")
    con <- dbConnect(drv, dbname=db, host=url, port=port, user=usrpwd[1],
                     password=usrpwd[2])
    out <- dbGetQuery(con, query)
    dbDisconnect(con)
    out
}
