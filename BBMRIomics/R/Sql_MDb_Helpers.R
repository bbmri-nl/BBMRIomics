##' get a view from the metadatabase
##'
##' @title get a view from the metadatabase
##' @param view Name of the to be retrieved view
##' @param ... arguments to be passed to \link[BBMRIomics]{runQuery}
##' @importFrom RPostgreSQL PostgreSQL dbConnect dbGetQuery dbDisconnect
##' @export
##' @return A data.frame containing the view.
##' @author Davy Cats
##' @examples
##' \dontrun{ 
##' view <- getSQLview("minimalphenotypes")
##' }
getSQLview <- function(view, verbose=T, ...) {
    views <- listViews(verbose=F, ...)
    
    if (! tolower(view) %in% views){
        stop(view, " is not a known view!")
    }
    
    out <- runQuery(paste0("SELECT * FROM ", view), verbose=verbose, ...)
    out
}


##' retrieve a list of available views
##'
##' @title retrieve a list of available views
##' @param ... arguments to be passed to \link[BBMRIomics]{runQuery}
##' @author Davy Cats
##' @importFrom RPostgreSQL PostgreSQL dbConnect dbGetQuery dbDisconnect
##' @export
##' @return A vector containing the view names.
##' @examples
##' \dontrun{
##' tables <- listViews()
##' }
listViews <- function(...) {
    out <- runQuery(paste0("SELECT table_name ",
                            "FROM information_schema.views ",
                            "WHERE table_schema = 'views'"),
                    ...)
    out$table_name
}

##' retrieve a list of available tables
##'
##' @title retrieve a list of available tables
##' @param ... arguments to be passed to \link[BBMRIomics]{runQuery}
##' @author Davy Cats
##' @importFrom RPostgreSQL PostgreSQL dbConnect dbGetQuery dbDisconnect
##' @export
##' @return A vector containing the table names.
##' @examples
##' \dontrun{
##' tables <- listTables()
##' }
listTables <- function(...) {
    out <- runQuery(paste0("SELECT tablename ",
                            "FROM pg_catalog.pg_tables ",
                            "WHERE schemaname = 'tables'"),
                    ...)
    out$tablename
}

##' send a query to the database
##'
##' @title send a query to the database
##' @param query SQL query to be send to the database
##' @param usrpwd Username and password concatenated by a colon,
##' defaults to "guest:guest".
##' @param url URL through which the database can be accessed, 
##' defaults to "localhost".
##' @param port port to be used to connect to the database, defaults
##' to 5432.
##' @param db name of the database to be connected to, defaults to
##' "rp3_rp4_meta".
##' @author Davy Cats
##' @importFrom RPostgreSQL PostgreSQL dbConnect dbGetQuery dbDisconnect
##' @export
##' @return A data.frame with the query results.
##' @examples 
##' \dontrun{
##' visits <- runQuery("SELECT * FROM visit;", RP3_MDB_USRPWD)
##' }
runQuery <- function(query, usrpwd="guest:guest", url="localhost", port=5432,
                      db="rp3_rp4_meta", verbose=T){
    if (usrpwd == "anonymous") {
        if (verbose) {
            message("Accessing the database as 'guest' user, which has ",
                    "limited read access.")
        }
        usrpwd <- "guest:guest"
    }

    usrpwd <- strsplit(usrpwd, ":")[[1]]
    if (verbose) {
        message("Accessing the '", db, "' database at '", url, ":", 
                    port, "' as user '", usrpwd[1], "'.")
    }
    
    drv <- PostgreSQL()
    con <- dbConnect(drv, dbname=db, host=url, port=port, user=usrpwd[1],
                     password=usrpwd[2])
    out <- dbGetQuery(con, query)
    dbDisconnect(con)
    out
}


##' retrieve database version
##'
##' @title retrieve database version
##' @param db name of the database to be connected to, defaults to
##' "rp3_rp4_meta".
##' @param ... arguments to be passed to \link[BBMRIomics]{runQuery}
##' @author Davy Cats
##' @importFrom RPostgreSQL PostgreSQL dbConnect dbGetQuery dbDisconnect
##' @export
##' @return The commit hash for the git repository from which the database 
##' was built.
##' @examples
##' \dontrun{
##' v <- mdbVersion()
##' }
mdbVersion <- function(db="rp3_rp4_meta", ...){
    out <- runQuery(paste0("SELECT description ",
                           "FROM pg_shdescription ",
                           "JOIN pg_database ON objoid = pg_database.oid ",
                           "WHERE datname = '", db, "'"), db=db, ...)
    out$description
}
