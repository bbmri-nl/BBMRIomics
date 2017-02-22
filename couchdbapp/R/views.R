
view <- function(viewname, db="https://metadatabase.bbmrirp3-lumc.surf-hosted.nl:6984/bios-test/", converter="json2csv", usrpwd=NULL) {

    require(httr)
    require(tibble)
    require(readr)

    set_config(config(ssl_verifypeer = 0L, ssl_verifyhost = 0L))
    url <- file.path(db, "_design/couchdbapp/_list", converter, viewname)      
    usrpwd <- unlist(strsplit(usrpwd, ":"))
    response <- GET(url, authenticate(usrpwd[1], usrpwd[2]), verbose())
    ##parsing done automatically including dates
    ##also problemes are reported
    as.data.frame(content(response, type="text/csv", as = "parsed", na = "null"))
}

phenotypes <- view("phenotypes?key=%22RS%22", usrpwd="mvaniterson:m0l3p1")

dim(phenotypes)
head(phenotypes)
str(phenotypes)
