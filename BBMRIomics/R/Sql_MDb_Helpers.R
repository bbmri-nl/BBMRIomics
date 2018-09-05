##' get a view from the metadatabase
##'
##' @title get a view from the metadatabase
##' @param view The name of the to be retrieved view.
##' @param usrpwd The username and password concatenated by a colon.
##' This defaults to "guest:guest".
##' @param url The URL through which the database can be accessed. This
##' defaults to "localhost".
##' @param port The port to be used to connect to the database. This defaults
##' to 5432.
##' @param db The name of database to be connected to. This defaults to
##' "rp3_rp4_meta".
##' @importFrom RPostgreSQL PostgreSQL dbConnect dbGetQuery dbDisconnect
##' @export
##' @return A data.frame containing the view.
##' @author Davy Cats
##' @examples
##' \dontrun{ 
##' view <- getSQLview("minimalphenotypes")
##' }
getSQLview <- function(view, usrpwd="guest:guest", url="localhost", port=5432,
                       db="rp3_rp4_meta") {
    views <- getViews(usrpwd=usrpwd, url=url, port=port, db=db)
    
    if (! tolower(view) %in% views){
        stop(view, " is not a known view!")
    }
    
    usrpwd <- strsplit(usrpwd, ":")[[1]]
    drv <- PostgreSQL()
    con <- dbConnect(drv, dbname=db, host=url, port=port, user=usrpwd[1],
                     password=usrpwd[2])
    tbl <- dbGetQuery(con, paste0("SELECT * FROM ", view))
    dbDisconnect(con)
    tbl
}


##' retrieve view names
##'
##' @title retrieve view names
##' @param usrpwd The username and password concatenated by a colon.
##' This defaults to "guest:guest".
##' @param url The URL through which the database can be accessed. This
##' defaults to "localhost".
##' @param port The port to be used to connect to the database. This defaults
##' to 5432.
##' @param db The name of database to be connected to. This defaults to
##' "rp3_rp4_meta".
##' @author Davy Cats
##' @importFrom RPostgreSQL PostgreSQL dbConnect dbGetQuery dbDisconnect
##' @export
##' @return A vector containing the view names.
##' @examples
##' \dontrun{
##' tables <- getViews()
##' }
getViews <- function(usrpwd="guest:guest", url="localhost", port=5432,
                     db="rp3_rp4_meta") {
    out <- .runQuery(paste0("SELECT table_name ",
                            "FROM information_schema.views ",
                            "WHERE table_schema = 'views'"),
                     usrpwd=usrpwd, url=url, port=port, db=db)
    out$table_name
}

##' retrieve table names
##'
##' @title retrieve table names
##' @param usrpwd The username and password concatenated by a colon.
##' This defaults to "guest:guest".
##' @param url The URL through which the database can be accessed. This
##' defaults to "localhost".
##' @param port The port to be used to connect to the database. This defaults
##' to 5432.
##' @param db The name of database to be connected to. This defaults to
##' "rp3_rp4_meta".
##' @author Davy Cats
##' @importFrom RPostgreSQL PostgreSQL dbConnect dbGetQuery dbDisconnect
##' @export
##' @return A vector containing the table names.
##' @examples
##' \dontrun{
##' tables <- getTables()
##' }
getTables <- function(usrpwd="guest:guest", url="localhost", port=5432,
                     db="rp3_rp4_meta") {
    out <- .runQuery(paste0("SELECT tablename ",
                            "FROM pg_catalog.pg_tables ",
                            "WHERE schemaname = 'tables'"),
                     usrpwd=usrpwd, url=url, port=port, db=db)
    out$tablename
}

##' send a query to the database
##'
##' @title send a query to the database
##' @param query The query to be send to the database.
##' @param usrpwd The username and password concatenated by a colon.
##' This defaults to "guest:guest".
##' @param url The URL through which the database can be accessed. This
##' defaults to "localhost".
##' @param port The port to be used to connect to the database. This defaults
##' to 5432.
##' @param db The name of database to be connected to. This defaults to
##' "rp3_rp4_meta".
##' @author Davy Cats
##' @importFrom RPostgreSQL PostgreSQL dbConnect dbGetQuery dbDisconnect
##' @return A data.frame with the query results.
##' @examples 
##' \dontrun{
##' visits <- .runQuery("SELECT * FROM visit;", RP3_MDB_USRPWD, RP3_MDB)
##' }
.runQuery <- function(query, usrpwd="guest:guest", url="localhost", port=5432,
                      db="rp3_rp4_meta"){
    usrpwd <- strsplit(usrpwd, ":")[[1]]
    if (usrpwd[1] == "guest") {
        message("Accesing the database as 'guest' user, which was limited ",
            "read access.")   
    }
    
    drv <- PostgreSQL()
    con <- dbConnect(drsv, dbname=db, host=url, port=port, user=usrpwd[1],
                     password=usrpwd[2])
    out <- dbGetQuery(con, query)
    dbDisconnect(con)
    out
}
