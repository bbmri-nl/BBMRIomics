##RP3/RP4 basic metadata overviews
##Morris Swertz
##Fri Mar 17 09:02:46 2017"

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
##subjects <- subset(subjects, !is.na(bios_id))
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
colnames(pheno) <- c("bios_id", "biobank_abbrv", "sex", "age")

sex <- pheno$sex
pheno$sex[sex] <- "male"
pheno$sex[!sex] <- "female"

##Fix LLS name
pheno$biobank_abbrv[grepl("LLS_", pheno$biobank_abbrv)] <- "LLS"
##Fix RS name
pheno$biobank_abbrv[grepl("ERF_ERGO", pheno$biobank_abbrv)] <- "RS"
##Fix LL name
pheno$biobank_abbrv[grepl("LIFELINES", pheno$biobank_abbrv)] <- "LL"
##Fix NTRT name
pheno$biobank_abbrv[grepl("VUNTR", pheno$biobank_abbrv)] <- "NTR"

head(pheno)

table(pheno$biobank_abbrv)

write.table(pheno, file="RP4_basic_overview.csv", row.names=FALSE, quote=FALSE, sep=",")

head(view)
head(pheno)

merged <- merge(view, pheno, by="bios_id", suffixes=c(".rp3", ".rp4"))

write.table(merged, file="merged_rp3rp3.csv", row.names=FALSE, quote=FALSE, sep=",")

##fix these relations once mdb is running again!!!
##Fri Mar  3 13:02:43 2017
##M van Iterson
relations <- getView("getRelations", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB)
complrel <- relations[!is.na(relations$relation_type), ]
setdiff(complrel$ids, complrel$relation_id)
setdiff(complrel$relation_id, complrel$ids)


##linking GoNL phenotypes to fast-files
##Harmen Draisma
##Fri Mar 10 14:12:32 2017
library(BBMRIomics)
ids <- getView("getIds", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB)
gonl <- subset(ids, !is.na(gonl_id))
dim(gonl)
##[1] 694  14

##use harmonized phenotypes
source(file.path(path.package("BBMRIomics"), "scripts/Phenotype_Helpers.R"), verbose=FALSE)
phenotypes <- cleanPhenotypes()
gonl <- merge(gonl, phenotypes, by="ids", all.x=TRUE)
dim(gonl)
##[1] 694  70

fastq <- getView("getFastq", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB)
gonl <- merge(gonl, fastq, by="ids")
dim(gonl) 
##[1] 516  74

##not for all gonl_id's there is RNAseq data available
##some individuals may have multiple fastq-files!

head(gonl)

files <- c(gonl$R1[1], gonl$R2[1])
srm.path <- file.path(dirname(files[1])) 
srm.path <- gsub("srm.*nl", SRM_BASE, srm.path) ##for curl access with need slightly different url                      
files <- basename(files)
srm.path
files
SRM2VM(srm.path, files, vm.path="~", proxy=GRID_PROXY)

##Rick/Bas
##NTR logitudinal samples
##Thu Mar 16 09:03:37 2017
library(BBMRIomics)

ids <- getView("getIds", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB)
ntr <- subset(ids, biobank_id == "NTR")

tbl <- table(ntr$person_id)
reps <- subset(ntr, person_id %in% names(tbl[tbl == 2]))
reps[order(reps$bios_id),c(1,5,8,9)]

dim(reps)

rels <- getView("getRelations", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB)
rels <- subset(rels, relation_type == "has repeated measurements")
dim(rels)
rels[order(rels$ids),]

sort(setdiff(unique(reps$ids), rels$ids))

pheno <- getView("allPhenotypes", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB)

mid <- match(rels$ids, pheno$ids)
rels$DNA_BloodSampling_Age <- pheno$DNA_BloodSampling_Age[mid]

rels$person_id <- unlist(lapply(as.character(rels$ids), function(x) unlist(strsplit(x, "-"))[2]))

dups <- rels$person_id[duplicated(rels[,c("person_id", "DNA_BloodSampling_Age")])]

missing <- rels$person_id[is.na(rels$DNA_BloodSampling_Age)]

subset(rels, person_id %in% dups | person_id %in% missing)


rels[order(rels$person_id),]

write.table(rels, file="NTR_longitudinal_ambiguous.csv", row.names=FALSE, quote=FALSE, sep=",")
