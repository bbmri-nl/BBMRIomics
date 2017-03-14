##RP3/RP4 basic metadata overviews
##Morris Swertz
##Tue Feb 28 14:36:28 2017

library(BBMRIomics)

##RP3
url <- "https://metadatabase.bbmrirp3-lumc.surf-hosted.nl:6984/bios/_design/overview/_view/basic"
usrpwd <- RP3_MDB_USRPWD
request <- paste0("curl -X GET ", url, " -u ",  usrpwd, " -k -g")
response <- system(request, intern = TRUE)
view <- jsonlite::fromJSON(response)$rows$value
head(view)
dim(view)
summary(view)
table(view$sex, useNA="always")

write.table(view, file="RP3_basic_overview.csv", row.names=FALSE, quote=FALSE, sep=",")

##RP4
molgenis.connect(usrpwd=RP4_DB_USRPWD, url=RP4_DB)
subjects <- molgenis.get.all("subjects", verbose=FALSE)
dim(subjects)
head(subjects)

##adding the real bios id
subjects <- subset(subjects, !is.na(bios_id))
table(subjects$biobank)

subjects$real_bios_id <- ""
ll <- subjects$biobank == "LIFELINES"
subjects$real_bios_id[ll] <- subjects$bios_id[ll] ##everything fine
rs <- subjects$biobank == "ERF_ERGO"
subjects$real_bios_id[rs] <- subjects$bios_id[rs] ##everything fine
lls <- subjects$biobank == "LLS_PARTOFFS" | subjects$biobank ==  "LLS_SIBS"
subjects$real_bios_id[lls] <- paste0("LLS-", subjects$bios_id[lls]) ##a little bit pasting
ids <- getView("getIds", url=RP3_MDB, usrpwd=RP3_MDB_USRPWD)
ntr <- subjects$biobank == "VUNTR"
id <- match(subjects$bios_id[ntr], ids$person_id[ids$biobank_id == "NTR"]) ##not bios id but person id
subjects$real_bios_id[ntr] <- ids$bios_id[ids$biobank_id == "NTR"][id]
subjects$real_bios_id[is.na(subjects$bios_id)] <- NA

table(subjects$biobank[!is.na(subjects$bios_id)])

samples <- molgenis.get.all("samples", verbose=FALSE)
colnames(subjects)[2] <- "subject_id" ##I find the naming a bit confusing
colnames(samples)[3] <- "sample_id"
id <- match(subjects$subject_id, samples$subject_id)
pheno <- cbind(subjects, samples[id,])

pheno <- pheno[, c("real_bios_id", "biobank", "gender", "age_at_collection")]
colnames(pheno) <- c("bios_id", "biobank_abbr_rp4", "sex", "age")
pheno$biobank_abbr_rp3 <- gsub("-.*$", "", pheno$bios_id)

sex <- pheno$sex
pheno$sex[sex] <- "male"
pheno$sex[!sex] <- "female"

head(pheno)

write.table(pheno, file="RP4_basic_overview.csv", row.names=FALSE, quote=FALSE, sep=",")

head(view)
head(pheno)

merged <- merge(view, pheno, by="bios_id", suffixes=c(".rp3", ".rp4"))

write.table(merged, file="merged_rp3rp3.csv", row.names=FALSE, quote=FALSE, sep=",")
