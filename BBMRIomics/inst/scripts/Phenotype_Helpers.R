##validate/harmonize phenotypes

.time_conversion <- function(x) {
    suppressPackageStartupMessages(require(lubridate))   
    suppressWarnings({
        y <- hms(x)
        y[is.na(y)] <- hm(x[is.na(y)])
    })
    as.character(y)
}

.date_conversion <- function(x) {
    suppressPackageStartupMessages(require(lubridate))
    suppressWarnings({
        y <- ymd(x)
        y[is.na(y)] <- dmy(x[is.na(y)])
    })
    as.character(y)
}

.rename <- function(phenotypes) {

    phenotypes$DNA_Extraction_Method <- gsub("salting out|Saltingout", "Salting out", phenotypes$DNA_Extraction_Method)
    phenotypes$DNA_Extraction_Method <- gsub("QIAamp DNA minikit", "DNA Mini Kit (Qiaamp)", phenotypes$DNA_Extraction_Method)

    phenotypes$RNA_Extraction_Method[!is.na(phenotypes$RNA_Extraction_Method)] <-"PAXgene Blood RNA Kit (Qiagen)"

    phenotypes$RNA_Source[!is.na(phenotypes$RNA_Source)] <- "Whole peripheral blood/PAX gene"
    phenotypes$RNA_Source[grepl("NTR", phenotypes$ids)] <- "Whole peripheral blood/PAX gene/Transferred from EDTA tube"

    phenotypes$GWAS_Chip <- gsub(" *, *", "/", phenotypes$GWAS_Chip)
    phenotypes$GWAS_Chip <- gsub("Illumina human omni express", "Illumina OmniExpress", phenotypes$GWAS_Chip)
    phenotypes$GWAS_Chip <- gsub("Illumina660", "Illumina Human660-Quad", phenotypes$GWAS_Chip)
    phenotypes$GWAS_Chip <- gsub("OverlappingSNPsfrom", "Overlapping SNPs from ", phenotypes$GWAS_Chip)
    phenotypes$GWAS_Chip <- gsub("GONLSequenceData", "GONL", phenotypes$GWAS_Chip)

    phenotypes$DNA_QuantificationMethod <- gsub("nanodrop|nano", "Nanodrop", phenotypes$DNA_QuantificationMethod)
    phenotypes$DNA_QuantificationMethod <- gsub("spectrofotometer|Spectofotometer|spectofotometer", "Spectrofotometer", phenotypes$DNA_QuantificationMethod)

    phenotypes$Ascertainment_criterion <- gsub("^GoNL$", "GONL_subject", phenotypes$Ascertainment_criterion)
    phenotypes$Ascertainment_criterion <- gsub("^Complete genomics sequencing$", "CG_subject", phenotypes$Ascertainment_criterion)
    phenotypes$Ascertainment_criterion <- gsub("^GoNL / Complete genomics sequencing", "GoNL/CG_subject", phenotypes$Ascertainment_criterion)

    phenotypes
}

.relabel <- function(phenotypes) {

    phenotypes$Sex <- factor(phenotypes$Sex)
    levels(phenotypes$Sex) <- c("Male", "Female")
    phenotypes$Sex <- as.character(phenotypes$Sex)

    phenotypes$Smoking <- factor(phenotypes$Smoking)
    levels(phenotypes$Smoking) <- c("non-smoker", "former smoker", "current smoker")
    phenotypes$Smoking <- as.character(phenotypes$Smoking)

    phenotypes$LipidMed <- factor(phenotypes$LipidMed)
    levels(phenotypes$LipidMed) <- c("no", "statins", "yes but no statins")
    phenotypes$LipidMed <- as.character(phenotypes$LipidMed)

    phenotypes$LDLcholMethod <- factor(phenotypes$LDLcholMethod)
    levels(phenotypes$LDLcholMethod) <- "Friedewald estimation"
    phenotypes$LDLcholMethod <- as.character(phenotypes$LDLcholMethod)

    phenotypes$Lipids_BloodSampling_Fasting <- factor(phenotypes$Lipids_BloodSampling_Fasting)
    levels(phenotypes$Lipids_BloodSampling_Fasting) <- c("no", "yes")
    phenotypes$Lipids_BloodSampling_Fasting <- as.character(phenotypes$Lipids_BloodSampling_Fasting)

    phenotypes
}

.addUnits <- function(phenotypes){
    path <- "/media/mvaniterson/Storage/packages/BIOSRutils/extdata"
    info <- read.table(file.path(path, "phenotypes.csv"), sep="\t", header=TRUE, as.is=TRUE)
    info <- info[!is.na(info$unit),]
    mid <- match(colnames(phenotypes), info$phenotype)
    colnames(phenotypes)[mid] <- paste0(colnames(phenotypes)[mid], " (", info$unit, ")")
    phenotypes
}


.reduceRelations <- function(colData) {

        relation_type <- tapply(colData$relation_type, colData$run_id, function(x) {
            if(all(is.na(x)))
                NA
            else
                paste0(x, collapse=",")
        })

        relation_uuid <- tapply(colData$relation_uuid, colData$run_id, function(x) {
            if(all(is.na(x)))
                NA
            else
                paste0(x, collapse=",")
        })

        colData <- colData[!duplicated(colData$run_id), ]
        colData$relation_uuid <- relation_uuid
        colData$relation_type <- relation_type
        colData
    }


cleanPhenotypes <- function() {
    require(BBMRIomics)
    
    message("Obtain phenotypes ... ")

    phenotypes <- getView("allPhenotypes", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB)

    message("Rename certain variables ... ")
    phenotypes <- .rename(phenotypes)

    message("Relabel certain variables ... ")
    phenotypes <- .relabel(phenotypes)

    message("Create single Sampling Age, Time and Date columns and harmonise ... ")

    ##create one sampling Age column
    Age <- phenotypes$DNA_BloodSampling_Age
    phenotypes <- phenotypes[, !grepl("Sampling_Age", colnames(phenotypes))]
    phenotypes$Sampling_Age <- Age

    ##Harmonize and create one sampling Time column
    phenotypes[, grep("Sampling_Time", colnames(phenotypes))] <- apply(phenotypes[, grep("Sampling_Time", colnames(phenotypes))], 2, .time_conversion)
    phenotypes[, grep("Sampling_Time", colnames(phenotypes))] <- apply(phenotypes[, grep("Sampling_Time", colnames(phenotypes))], 2, function(x) gsub("^8H -11M 0S$", "9H 30M 0S", x)) ##the average
    phenotypes[, grep("Sampling_Time", colnames(phenotypes))] <- apply(phenotypes[, grep("Sampling_Time", colnames(phenotypes))], 2, function(x) gsub("^0S$", NA, x))

    Time <- phenotypes$DNA_BloodSampling_Time
    phenotypes <- phenotypes[, !grepl("Sampling_Time", colnames(phenotypes))]
    phenotypes$Sampling_Time <- Time

    ##Harmonize and create one sampling Date column
    phenotypes[, grep("Sampling_Date", colnames(phenotypes))] <- apply(phenotypes[, grep("Sampling_Date", colnames(phenotypes))], 2, .date_conversion)

    Date <-  phenotypes$DNA_BloodSampling_Date
    phenotypes <- phenotypes[, !grepl("Sampling_Date", colnames(phenotypes))]
    phenotypes$Sampling_Date <- Date

    missing <- 100*apply(phenotypes, 2, function(x) sum(is.na(x)))/nrow(phenotypes)

    message("Drop columns completely empty ...", paste0(colnames(missing[missing > 95]), collapse=", "))

    phenotypes <- phenotypes[, !(colnames(phenotypes) %in% c("CHCM", "CH","HDW"))]

    ##Calculate LDLchol, for those missing, from Tot, HDL and TriGlycerides
    ##http://www.gpnotebook.co.uk/simplepage.cfm?ID=x20030114211535665170
    LDLchol <- phenotypes$TotChol - phenotypes$HDLchol - (phenotypes$Triglycerides/2.2)
    phenotypes$LDLchol[is.na(phenotypes$LDLchol)] <- LDLchol[is.na(phenotypes$LDLchol)]

    ##phenotypes <- .addUnits(phenotypes)
    
    invisible(phenotypes)
}




if(FALSE) {

    .is.failed <- function(validatorObject){
        apply(values(validatorObject), 1, function(x) {
            x[is.na(x)] <- TRUE
            !all(x)
        })
    }


    ##--------------------
    ##Validation
    ##-------------------

    library(BBMRIomics)
    library(validate)

    phenotypes <- getView("allPhenotypes", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB)

    ##validate categorical variables
    cat <- check_that(phenotypes,
                      Sex == 0 | Sex == 1,
                      Smoking  == 0 | Smoking == 1 | Smoking == 2,
                      LipidMed == 0 | LipidMed == 1 | LipidMed == 2,
                      LDLcholMethod == 1 | LDLcholMethod == 2,
                      Lipids_BloodSampling_Fasting == 0 | Lipids_BloodSampling_Fasting == 1)

    summary(cat)

    ##validate continuous variables
    con <- check_that(phenotypes,
                      Height > 0 & Height < 250,
                      Weight > 0 & Weight < 150,
                      TotChol > 0,
                      LDLchol > 0,
                      HDLchol > 0,
                      Triglycerides > 0,
                      hsCRP > 0,
                      LUC > 0,
                      Mono > 0,
                      Eos >= 0, ##really???
                      Lymph > 0,
                      Neut > 0,
                      Baso >= 0,
                      Granulocyte > 0,
                      CHCM > 0,
                      WBC > 0,
                      RBC > 0,
                      RDW > 0,
                      CH > 0,
                      HCT > 0,
                      HGB > 0,
                      PLT > 0,
                      HDW > 0,
                      MCHC > 0,
                      MPV > 0,
                      MCH > 0,
                      MCV > 0,
                      RNA_RIN > 0,
                      BirthYear > 1900 & BirthYear < 2016)

    summary(con)

    phenotypes[is.failed(con), c("Eos", "Baso")]

    ##validate ratio/percentages
    rp <- check_that(phenotypes,
                     RNA_A260280ratio != 0,
                     DNA_A260A280ratio != 0,
                     Granulocyte_Perc > 0 & Granulocyte_Perc < 100,
                     LUC_Perc > 0 &  LUC_Perc < 100,
                     Mono_Perc > 0 & Mono_Perc < 100,
                     Eos_Perc >= 0 & Eos_Perc < 100,
                     Lymph_Perc > 0 & Lymph_Perc < 100,
                     Neut_Perc > 0 & Neut_Perc < 100,
                     Baso_Perc >= 0 & Baso_Perc < 100,
                     Mono_Perc + Lymph_Perc + Neut_Perc + Baso_Perc +  Eos_Perc > 90)

    summary(rp)

    ##validate columns containing Age variable
    grep("Age", colnames(phenotypes), value=TRUE)

    age <- check_that(phenotypes,
                      var_group(DNA_BloodSampling_Age, Lipids_BloodSampling_Age,
                                RNA_BloodSampling_Age, LipidsMed_Age,
                                CRP_BloodSampling_Age, Anthropometry_Age,
                                CellCount_BloodSampling_Age, Smoking_Age) > 0,
                      var_group(DNA_BloodSampling_Age, Lipids_BloodSampling_Age,
                                RNA_BloodSampling_Age, LipidsMed_Age,
                                CRP_BloodSampling_Age, Anthropometry_Age,
                                CellCount_BloodSampling_Age, Smoking_Age) < 125,
                      rowSds(cbind(DNA_BloodSampling_Age, Lipids_BloodSampling_Age, RNA_BloodSampling_Age,
                                   LipidsMed_Age, CRP_BloodSampling_Age, Anthropometry_Age, CellCount_BloodSampling_Age, Smoking_Age), na.rm=TRUE) == 0)

    summary(age)

    table(gsub("-.*$", "", phenotypes$ids)[!is.failed(age) & !is.na(phenotypes$Smoking_Age)])

    pairs(phenotypes[is.failed(age), grep("Age", colnames(phenotypes))], pch=".")

    ##harmonize times
    phenotypes[, grep("Time", colnames(phenotypes))] <- apply(phenotypes[, grep("Time", colnames(phenotypes))], 2, time_conversion)

    phenotypes[, grep("Time", colnames(phenotypes))] <- apply(phenotypes[, grep("Time", colnames(phenotypes))], 2, function(x) gsub("^8H -11M 0S$", "9H 30M 0S", x)) ##the average

    phenotypes[, grep("Time", colnames(phenotypes))] <- apply(phenotypes[, grep("Time", colnames(phenotypes))], 2, function(x) gsub("^0S$", NA, x))

    grep("Time", colnames(phenotypes), value=TRUE)

    times <- rowSums(with(phenotypes, cbind(CellCount_BloodSampling_Time,
                                            CRP_BloodSampling_Time,
                                            RNA_Sampling_Time,
                                            Lipids_BloodSampling_Time) == DNA_BloodSampling_Time), na.rm=TRUE)

    table(times)
    phenotypes[times==0, grep("Time", colnames(phenotypes))]
    phenotypes[times==3, grep("Time", colnames(phenotypes))]

    ##harmonize dates
    grep("Sampling_Date", colnames(phenotypes), value=TRUE)

    phenotypes[, grep("Sampling_Date", colnames(phenotypes))] <- apply(phenotypes[, grep("Sampling_Date", colnames(phenotypes))], 2, date_conversion)

    head(phenotypes[, grep("Sampling_Date", colnames(phenotypes))])

    dates <- apply(phenotypes[, grep("Sampling_Date", colnames(phenotypes))], 1, function(x)  sum(duplicated(x[!is.na(x)])))

    table(dates)
    phenotypes[dates==0, grep("Sampling_Date", colnames(phenotypes))]
    phenotypes[dates==1, grep("Sampling_Date", colnames(phenotypes))]
    phenotypes[dates==3, grep("Sampling_Date", colnames(phenotypes))]
    phenotypes[dates==4, grep("Sampling_Date", colnames(phenotypes))]

}
