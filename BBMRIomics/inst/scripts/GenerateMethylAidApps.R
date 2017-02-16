##This script generates the shiny MethylAid apps

library(BIOSRutils)
APPDIR <- "/srv/shiny-server"

for(biobank in BIOBANKS) {

    datafile <- dir(file.path(RP3DATADIR, "IlluminaHumanMethylation450k", "450k", biobank), full.names=TRUE, pattern="MethylAid.RData")

    message(datafile)    

    dir.create(file.path(APPDIR, "BIOSRutils/DNAm", biobank), recursive=TRUE)
    
    cat("require(shiny)\n",
        "require(MethylAid)\n",
        "require(MethylAidData)\n",
        "data(exampleDataLarge)\n",
        paste0("load('", datafile, "')"), "\n",
        "shinyApp(ui = MethylAid:::ui450k(MethylAid),
                  server = MethylAid:::server450k(MethylAid,
                                         thresholds = list(MU = 10.50, OP = 11.75, BS = 12.75, HC = 13.25, DP = 0.95), background=exampleDataLarge))",
         file=file.path(APPDIR, "BIOSRutils/DNAm", biobank, "App.R"))
        
    }


