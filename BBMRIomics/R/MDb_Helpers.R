##' get view from metadatabase
##'
##' extract information from the metadatabase using predefined views
##' and a list function to convert json to csv.xx
##' IN DEVELOPMENT using list functions and the readr-package
##' @title get view from metadatabase
##' @param viewname name of the view
##' @param db metadatabase url
##' @param converter default json2csv
##' @param usrpwd metadata usrpwd
##' @return metadatabase view converted to a data.frame
##' @author mvaniterson
##' @importFrom readr read_csv
##' @importFrom httr set_config config GET authenticate content
view <- function(viewname, db="https://metadatabase.bbmrirp3-lumc.surf-hosted.nl:6984/bios-test/", converter="json2csv", usrpwd=NULL) {

    set_config(config(ssl_verifypeer = 0L, ssl_verifyhost = 0L))
    url <- file.path(db, "_design/couchdbapp/_list", converter, viewname)
    request <- paste0(url, "?reduce=false")
    usrpwd <- unlist(strsplit(usrpwd, ":"))
    response <- GET(url, authenticate(usrpwd[1], usrpwd[2])) ##verbose()
    ##parsing done automatically including dates
    ##also problemes are reported
    as.data.frame(content(response, type="text/csv", as = "parsed", na="null"))
}

##' get view from metadatabase
##'
##' TODO handle reduce options properly
##' @title get view from metadatabase
##' @param view name of the view
##' @param url unique resource location of the metadatabase
##' @param usrpwd username and passwrd concatenated by colon 'usr:pwd'
##' defaults to 'anonymous' account will
##' be used using stored output of the metadatabase
##' @param selection some view have a reduce option
##' @param force.keys FALSE
##' @param row.names use ids as row names in the returned data.frame
##' @param verbose default TRUE prints some additional output
##' @importFrom jsonlite fromJSON toJSON
##' @importFrom utils read.table
##' @export
##' @return data.frame containing the view
##' @author mvaniterson
##' @examples
##' \dontrun{ 
##' view <- getView("getIds", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB)
##' head(view)
##' colnames(view)
##' dim(view)
##' view <- getView("freeze1RNASeq", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB)
##' head(view)
##' colnames(view)
##' dim(view)
##' view <- getView("allPhenotypes", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB)
##' head(view)
##' dim(view)
##' ##For example make a plot of the Age distribution in BIOS
##' view <- getView("minimalPhenotypes", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB)
##' library(ggplot2)
##' view <- view[, c("biobank_id", "sex", "DNA_BloodSampling_Age")]
##' ##convert to numeric
##' view$DNA_BloodSampling_Age <- as.numeric(view$DNA_BloodSampling_Age)
##' str(as.numeric(view$DNA_BloodSampling_Age))
##' g <- ggplot(view, aes(x = DNA_BloodSampling_Age))
##' g <- g + geom_histogram()
##' g <- g + facet_grid(. ~ biobank_id)
##' g
##' g <- g + facet_grid(sex ~ biobank_id)
##' g
##' }
getView <- function(view, url, usrpwd="anonymous", selection="?reduce=false", force.keys=FALSE, row.names=FALSE, verbose = FALSE){
    
    ##views/designs
    views <- list()
    views[["Metadatabase"]] <- list(EGA = c("freeze1RNASeq","freeze1Methylation", "freeze2RNASeq","freeze2Methylation"),
                  Files = c("getFastq", "getIdat"),
                  Identifiers = c("getIds"),
                  Phenotypes = c("allPhenotypes", "cellCounts", "minimalPhenotypes"),
                  Runs = c("getGenotypes", "getMethylationRuns",
                           "getRNASeqRuns", "getImputations",
                           "getRelations"),
                  Samplesheets = c("rnaseqSamplesheet", "methylationSamplesheet"),
                  Verification = c("md5sum"))
    
    views[["Rundatabase"]] <- list(stats=c("getStats", "getBAM"))
    
    getDesign <- function(view) {
        ret <- names(which(rapply(views, function(x) x == view)))
        dsgn <- gsub("^.*\\.|[0-9]$", "", ret)
        names(dsgn) <- gsub("\\..*$", "", ret)
        dsgn
    }

    design <- getDesign(view)
    
    if(length(design) != 1)
        stop(view, "not known!")
 
    if(usrpwd != "anonymous") {
        view <- paste("_design", design , "_view", view, sep="/")
        if(!is.null(selection))
            view <- paste0(view, selection)

        ##option -g url encoding
        ##somethimes -3 is necessary for ssl3 protocol
        request <- paste0("curl -X GET ", url, view, " -u ", usrpwd, " -k -g")
        if(verbose)
            message(request)

        response <- system(request, intern=TRUE)
        if(any(grepl("error|Error", response)))
            stop(response)
    }
    else {
        message("No username and password provided for the MDB use stored views!")
        response <- file.path(path.package("BBMRIomics"), "views", paste0(design, "_", view, ".json"))
    }

    view <- fromJSON(response)

    if(verbose & !is.null(view$total_rows))
        message("Got ", view$total_rows, " records from database...")

    values <- view$rows$value
    ids <- view$rows$id
    keys <- view$rows$key

    if(is.null(ids) | force.keys)
        ids <- keys

    if(row.names)  {
        if(any(duplicated(ids)))
            stop("duplicated ids!")
        rownames(values) <- ids
    }
    else
        values <- cbind(ids, values)

    return(values)
}


##' extract document from couchdb database
##'
##' R-wrapper around curl to extract document from couchdb database
##' convert JSON to list of lists
##' @title extract document from couchdb database
##' @param id doc._id
##' @param url couchdb url
##' @param usrpwd username:password
##' @param verbose default TRUE
##' @return return list of lists
##' @author mvaniterson
##' @importFrom jsonlite fromJSON
##' @examples
##' \dontrun{
##' doc <- .getDoc("LLS-2749", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB)
##' names(doc)
##' str(doc)
##' }
.getDoc <- function(id, url, usrpwd, verbose=TRUE) {
    request <- paste0("curl -X GET ", url, id, " -u ", usrpwd, " -k -g")

    if(verbose)
        message(request)

    response <- system(request, intern=TRUE)
    if(any(grepl("error|Error", response)))
        stop(response)

    if(length(response) > 1)
        response <- paste0(response, collapse="")

    fromJSON(response)
}

##' put document in couchdb database
##'
##' converts R list of list to JSON (number of digits = 12)
##' and validate according to JSON-schema
##' and if validation succeeds puts document in database
##' @title put document in couchdb database
##' @param doc list of lists
##' @param url couchdb url
##' @param usrpwd username:password
##' @param verbose default TRUE
##' @param ... path of database schema
##' @return nothing
##' @author mvaniterson
##' @importFrom jsonlite toJSON
##' @examples
##' \dontrun{
##' doc <- .getDoc("LLS-2749", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB)
##' ##modify or add something
##' SCHEMA <- file.path(VM_BASE_ANALYSIS, "BBMRIomics/couchdbapp/schema/bios.json")
##' .validateDoc(doc, SCHEMA=SCHEMA)
##' ##BE REALLY CAREFUL WHAT YOU DO!!!
##' ##.putDoc(doc, usrpwd=RP3_MDB_USRPWD, url=RP3_MDB, SCHEMA=SCHEMA) 
##' }
.putDoc <- function(doc,  url, usrpwd, verbose=TRUE, ...) {

    json <- toJSON(doc, digits = 12, auto_unbox = TRUE) ##set number of digits large than default of 6
    json <- gsub("\\{\\}", "null", json) ## empty list set to "null"

    if(!is.null(.validateDoc(json, ...)))
        stop("Validation error!:", .validateDoc(json, ...))
    else
        message("Valid json document!")

    request <- paste0("curl -v -H 'Content-Type: application/json' -d ", shQuote(json), " -X PUT ", url, doc$`_id`, " -u ", usrpwd, " -k")

    if(verbose)
        message(request)

    system(request, intern=TRUE)
}

##' validate json document against schema
##'
##' uses jsonschema CLI Python script from Julian Berman:
##' https://github.com/Julian/jsonschema
##' @title validate json documents against schema
##' @param json R list representation of json document or json document
##' @param SCHEMA path to json schema: bios-schema/bios.json
##' @return error or NULL
##' @author mvaniterson
##' @examples
##' \dontrun{ 
##' ids <- getView("getIds", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB)
##' for(id in ids$ids) {
##'     print(id)
##'     doc <- .getDoc(id, usrpwd=RP3_MDB_USRPWD, url=RP3_MDB)
##'     val <- .validateDoc(doc, SCHEMA=file.path(VM_BASE_ANALYSIS, "BBMRIomics/couchdbapp/schema/bios.json"))
##'     if(!is.null(val))
##'         stop(val)
##' }
##' }
.validateDoc <- function(json, SCHEMA) {

    if(system2("jsonschema") == 127)
        stop("jsonschema probably not installed use: pip install git+https://github.com/Julian/jsonschema.git")
        
    if(class(json) == "list") {
        json <- toJSON(json, digits = 12, auto_unbox = TRUE, pretty=TRUE, na="null") ##set number of digits large than default of 6
        json <- gsub("\\{\\}", "null", json) ## empty list set to "null"
    }

    if(class(json) != "json")
        stop("json-object not of class 'json'!")

    doc.json <- tempfile(fileext=".json")

    doc.json.con <- file(doc.json)

    writeLines(json, con=doc.json.con)

    wd <- getwd()
    setwd(dirname(SCHEMA))
    on.exit(setwd(wd))

    request <- paste("jsonschema -V Draft4Validator -i", doc.json, SCHEMA)
    returned <- system(request, intern=TRUE)

    close(doc.json.con)
    return(attributes(returned)$status)
}
