##RP3/RP4 basic metadata overviews
##Morris Swertz
##Tue Mar 28 09:36:50 2017

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

colnames(view)[1] <- "rp3_id"
colnames(view)[5] <- "rp2_id"

##write.table(view, file="RP3_basic_overview.csv", row.names=FALSE, quote=FALSE, sep=",")

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

pheno <- pheno[, c("real_bios_id", "biobank", "gender", "age_at_collection", "subject_id")]
colnames(pheno) <- c("rp3_id", "biobank_abbrv", "sex", "age", "rp4_id")

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
pheno$metabolite <- TRUE

table(pheno$biobank_abbrv)

##write.table(pheno, file="RP4_basic_overview.csv", row.names=FALSE, quote=FALSE, sep=",")

view[1:10,]

pheno[1:10,]

merged <- merge(view, pheno, by=c("rp3_id", "biobank_abbrv"), all=TRUE, suffixes=c("_rp3", "_rp4"))

merged$metabolite[is.na(merged$metabolite)] <- FALSE
merged[1:10,]
dim(merged)

summary(merged)

age <- merged$age_rp3
age[is.na(age)] <- merged$age_rp4[is.na(age)]

sex <- merged$sex_rp4
sex[is.na(sex)] <- merged$sex_rp3[is.na(sex)]

merged$sex <- sex
merged$age <- age

merged$GoNL <- !is.na(merged$rp2_id)

merged <- merged[order(merged$biobank_abbrv),]
merged$uuid <- paste0("BBMRI_", 1:nrow(merged))

order <- c("uuid", "rp2_id", "rp3_id", "rp4_id", "biobank_abbrv", "GoNL", "DNA", "DNAm", "rnaseq", "metabolite", "sex", "age", "smoking", "wbcc")
merged <- merged[, order]

merged[1:10,]

merged$GoNL[is.na(merged$GoNL)] <- FALSE
merged$DNA[is.na(merged$DNA)] <- FALSE
merged$DNAm[is.na(merged$DNAm)] <- FALSE
merged$rnaseq[is.na(merged$rnaseq)] <- FALSE
merged$sex[is.na(merged$sex)] <- FALSE
merged$smoking[is.na(merged$smoking)] <- FALSE
merged$wbcc[is.na(merged$wbcc)] <- FALSE

merged$rp3_id <- NA
merged$rp4_id <- NA
merged[1:20,]

write.table(merged, file="merged_rp4rp3.csv", row.names=FALSE, quote=FALSE, sep=",")

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


##Joyce
##Update RS phenotypes
##Fri Apr 14 08:49:51 2017
library(BBMRIomics)
id <- getView("getIds", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB)
id <- subset(id, biobank_id == "RS")

rs <- getView("allPhenotypes", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB)
rs <- subset(rs, grepl("RS", ids))
rs <- droplevels(rs)



y <- sort(colnames(rs))

gonl <- read.table("16_12_2015_BBMRI_RP3_Rotterdam_rest_samples.csv", header=TRUE, sep="\t", na.strings=c(NA, NULL, "", "missing", "Missing"))
dim(gonl)
##fix columns names (manually)
x <- sort(colnames(gonl))
x <- x[-24]
mid <- match(x,y)
z <- cbind(x=x, y=y)
z[is.na(mid),]

##fix content
gonl$ids <- paste0("RS-", gonl$ids)
gonl$ids

gonl$Sex <- as.integer(gonl$Sex)
gonl$Sex[gonl$Sex == 2] <- 0
gonl$Sex

cname <- colnames(gonl)
i <- 6
cname[i]
gonl[,cname[i]]
rs[, colnames(rs) == cname[i]]


rs <- merge(rs, id, by="ids", all.y=TRUE)
rs$ids

mid <- match(gonl$ids, rs$ids)
rs[mid,]

##Adding Age's and Sex
ids <- gonl$ids

template <- BBMRIomics:::.getDoc("CODAM-2001", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB)

BBMRIomics:::.validateDoc(template, SCHEMA=file.path(VM_BASE_ANALYSIS, "BBMRIomics/couchdbapp/schema/bios.json"))

phenotypes <- names(template$phenotype)
phenotypes <- phenotypes[!grepl("Age|Sex", phenotypes)]

pheno <- as.list(rep(NA, length(phenotypes)))
names(pheno) <- phenotypes

length(pheno)

i <- 1

id <- ids[i]
message(id)
doc <- BBMRIomics:::.getDoc(id, usrpwd=RP3_MDB_USRPWD, url=RP3_MDB)
if(!is.null(doc$phenotype)) stop("has already phenotype info")

ph <- c(gonl[i, grepl("Age|Sex", colnames(gonl))], pheno)

doc$phenotype <- ph
doc$phenotype

BBMRIomics:::.validateDoc(doc, SCHEMA=file.path(VM_BASE_ANALYSIS, "BBMRIomics/couchdbapp/schema/bios.json"))


##Verify mismatches
##

rbind.all.columns <- function(x, y) {

    x.diff <- setdiff(colnames(x), colnames(y))
    y.diff <- setdiff(colnames(y), colnames(x))

    x[, c(as.character(y.diff))] <- NA

    y[, c(as.character(x.diff))] <- NA

    return(rbind(x, y))
}


path <- "/media/mvaniterson/Storage/beamer"

for(cohort in c("PAN", "CODAM", "RS", "LL", "LLS", "NTR")){

    files <- dir(file.path(path, "SwapDetection/"), pattern="ALL.txt$", full.names=TRUE)
    files <- c(files,  dir(file.path(path, "SwapDetection/"), pattern=paste0(cohort, "txt$"), full.names=TRUE))
    message(files)

    mm <- read.table(files[1], header=TRUE, sep="\t")
    mm$file <- basename(files[1])
    for(i in 2:length(files)) {
        m <- read.table(files[i], header=TRUE, sep="\t")
        m$file <- basename(files[i])
        mm <- rbind.all.columns(mm, m)
    }

    ptrn <- paste0("^",cohort,"-")
    view <- mm[grepl(ptrn, mm$ids.x) | grepl(ptrn, mm$ids.y) | grepl(ptrn, mm$colnames.x) | grepl(ptrn, mm$colnames.y), ]

    write.table(view[order(view$ids.x), ], file=file.path(path, paste0("SwapDetection/", cohort, "_view.txt")), sep=",", col.names=TRUE, row.names=FALSE)
}


##Remove invalid idat
id <- "9340996015_R06C02"
library(BBMRIomics)
runs <- getView("getMethylationRuns", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB)
subset(runs, grepl(id, run_id))


##something wrong with this one AC43B7ACXX-5-6?
runs <- getView("getRNASeqRuns", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB)

subset(runs, run_id == "AC43B7ACXX-5-6")

NTR-A969D-NT0027660-10392
bam <- getView("getBAM", usrpwd=RP3_MDB_USRPWD, url=RP3_RDB)

subset(bam, ids=="AC43B7ACXX-5-6")

bam$path[159]

file <- basename(bam$path[159])
srm.path <- file.path(dirname(bam$path[159]))
srm.path <- gsub("srm.*nl", SRM_BASE, srm.path) ##for curl access with need slightly different url
srm.path
file
SRM2VM(srm.path, file, vm.path="~", proxy=GRID_PROXY)


##CODAM genotype and sex
##Joost Verlouw
##Thu May 18 10:52:15 2017
library(BBMRIomics)
geno <- getView("getImputations", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB)
pheno <- getView("allPhenotypes", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB)

geno <- subset(geno, biobank_id == "CODAM")

geno2pheno <- merge(geno, pheno, all.x=TRUE, by="ids")
str(geno2pheno)

geno2pheno[,c("imputation_id", "imputation_reference", "Sex", "ids")]


##Check sex
##PB
library(BBMRIomics)

load(file.path(VM_BASE_DATA, "RNASeq/v2.1.3/gene_read/", paste0("rnaSeqData_ReadCounts_GoNL.RData")))

cnts <- assays(counts)$data[rownames(counts) %in% c("ENSG00000229807", "ENSG00000129824"),]

cex <- rep(1, ncol(counts))
cex[counts$gonl_id %in% c("gonl-186c", "gonl-228c", "gonl-231c", "gonl-24c")] <- 3
table(cex)
sex <- ifelse(counts$Sex == "Female", 1, 2)

pdf("~/gonlxistandrps4y1.pdf")
plot(log2(1+t(cnts)), col=adjustcolor(sex, alpha.f = 0.3), pch=16, cex=cex, ylab="XIST (ENSG00000229807)", xlab=" RPS4Y1 (ENSG00000129824)", main="log2(1+raw counts)")
dev.off()



##Red Blood Counts
##Rene Luijk
##Fri Jun 16 09:54:44 2017
library(BBMRIomics)
data(rnaSeqData_ReadCounts_BIOS_cleaned)
counts
library(edgeR)
y <- DGEList(assays(counts)$data, remove.zeros=TRUE)
y <- calcNormFactors(y)
y <- cpm(y, log=TRUE)
library(irlba)
pc <- prcomp_irlba(t(y), n=25)
summary(pc)
plot(pc$x[,c(1,2)], col=ifelse(colData(counts)$Sex == "Female", 1, 2))

library(pheatmap)

cov <- as.data.frame(colData(counts))[, c("RBC", "MCH", "MCV", "MPV", "HGB", "RDW", "HCT", "Eos_Perc", "Lymph_Perc", "Mono_Perc", "Neut_Perc")]

cxy1 <- cor(pc$x[, 1:25], cov, use="complete.obs")
cxy2 <- cor(cov, use="complete.obs")

pdf("inspect_RBC.pdf", onefile=TRUE)
pheatmap(cxy1)
pheatmap(cxy2)
dev.off()


##Adding HRCv1.1
##Joost Verlouw
##Tue Jul  4 2017

##data.frame with between bios-ids and HRCv1.1 imputation ID's

##this is an small example for adding one document
##for example adding HRCv1.1 to LLS-2447

library(BBMRIomics)
##functions that interact directly with the metadatabase and can
##modify data are not exported, for savety, therefore you need to use
##`:::` the triple-colon operator (which is not known to regular R users)

##extract the document (this could be based on your mapping file/data.frame)
doc <- BBMRIomics:::.getDoc("LLS-2447", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB)

##we need to add to files, genotype, imputation

str(doc$files$genotype$imputation) ##this is a data.frame

##extract the current information
imputation <- doc$files$genotype$imputation

##add new information (possibly extracted from your data.frame)
imputation <- rbind(imputation,
                    data.frame(reference="HRCv1.1",
                               id = "2447_2447"))
##check it
imputation

##replace old with new
doc$files$genotype$imputation <- imputation

##check again
doc$files$genotype$imputation

##a dry run only validation
SCHEMA <- file.path(VM_BASE_ANALYSIS, "BBMRIomics/couchdbapp/schema/bios.json")
BBMRIomics:::.validateDoc(doc, SCHEMA=SCHEMA)

##now add to database also validation is performed
BBMRIomics:::.putDoc(doc, usrpwd=RP3_MDB_USRPWD, url=RP3_MDB, SCHEMA=SCHEMA)
##if its returns ok your document is updated

##NOW IT IS REALLY ADDED EITHER USE FUTON TO REMOVE THE EXAMPLE OR
##RERUN AND REMOVE THE ADDED SECTION

RS

for(i in 1:nrow(RS)){
    biosid <- RS[i,1]
    id <- RS[i,2]
    doc <- BBMRIomics:::.getDoc(biosid, usrpwd=RP3_MDB_USRPWD, url=RP3_MDB)
    ##extract the current information
    imputation <- doc$files$genotype$imputation

    ##add new information (possibly extracted from your data.frame)
    imputation <- rbind(imputation,
                        data.frame(reference="HRCv1.1",
                                   id = id))   
    ##replace old with new
    doc$files$genotype$imputation <- imputation
    BBMRIomics:::.putDoc(doc, usrpwd=RP3_MDB_USRPWD, url=RP3_MDB, SCHEMA=SCHEMA)    

}



library(BBMRIomics)
gonl <- c("gonl-186c", "gonl-228c", "gonl-231c", "gonl-24c")
ids <- getView("getIds", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB)

colnames(ids)

gonl <- subset(ids, !is.na(gonl_id) | !is.na(old_gonl_id))

write.table(gonl[, 11:12], file="gonl_old_new.csv", sep="\t", row.names=FALSE, quote=FALSE)



##Bas
##Freeze2 numbers
##2017-07-13
library(BBMRIomics)
packageVersion("BBMRIomics")
##[1] ‘2.0.1’

data(package="BBMRIomics")

##1. Aantal unrelated met 450k
data(methData_Betas_BIOS_F2_cleaned)
ncol(betas)
##[1] 4386

##2. Aantal unrelated met 450k & RNA-seq
data(methData_Betas_BIOS_F2_cleaned)
data(rnaSeqData_ReadCounts_BIOS_cleaned)
length(intersect(colData(betas)$uuid, colData(counts)$uuid))
##[1] 3400

##3. Aantal unrelated met 450k & SNPs
sum(!is.na(colData(betas)$imputation_id))
##[1] 4025

##4. Aantal unrelated met 450k & RNA-seq & SNPs
betas <- betas[,!is.na(betas$imputation_id)]
ncol(betas) ##same as 3.
length(intersect(colData(betas)$uuid, colData(counts)$uuid))
##3219

##5. Aantal unrelated met RNA-seq & SNPs
data(rnaSeqData_ReadCounts_BIOS_cleaned)
sum(!is.na(colData(counts)$imputation_id))
##[1] 3357

##6. Aantal MZ tweeling paren met 450k
runs <- getView("getMethylationRuns", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB)
runs <- subset(runs, qc == "passed" & biobank_id == "NTR")
relations <- getView("getRelations", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB)
mz <- subset(relations, relation_type == "has monozygotic twin")
[1] 1881    7
mz <- mz[mz$ids %in% runs$ids & mz$relation_id %in% runs$ids,]
##[1] 1719    7 ##twin pairs

##7. Aantal MZ tweeling paren met RNA-seq
runs <- getView("getRNASeqRuns", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB)
runs <- runs[!duplicated(runs$run_id),] ##drop overlapping freezes
runs <- subset(runs, qc == "passed") ##keep those passing qc
runs <- subset(runs, type != "replicate") ##drop replicated
##now we have still original, reruns and merged runs
merged <- subset(runs, type == "merged" & freeze == 2) ##select
original <- subset(runs, !(ids %in% merged$ids) & type == "original")
runs <- rbind(merged, original)

runs <- subset(runs, qc == "passed" & biobank_id == "NTR")
relations <- getView("getRelations", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB)
mz <- subset(relations, relation_type == "has monozygotic twin")
dim(mz)
[1] 1881    7
mz <- mz[mz$ids %in% runs$ids & mz$relation_id %in% runs$ids,]
dim(mz)
##[1] 1090    7 ##twin pairs

##8. Aantal CpGs voor analyse dat overblijft na QC 450k arrays.
nrow(betas)
##this depends on the number of samples as well since we use filter criteria that depend on the number of samples!
##[1] 481388 ##no SNPs or sex-chromosomes removed


##EGA upload
library(BBMRIomics)
data(package="BBMRIomics")

##EGA_upload_unrelated_DNAm
data(methData_Betas_BIOS_F2_cleaned)

idat <- getView("getIdat", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB)
head(idat)

idat$idat <- apply(idat[, c("Red", "Grn")], 1, paste, collapse=", ")

metadata <- as.data.frame(colData(betas)[, c("Sex", "Sampling_Age", "Baso_Perc", "Eos_Perc", "Neut_Perc", "LUC_Perc", "Lymph_Perc", "Mono_Perc", "Baso", "Eos", "Neut", "LUC", "Lymph", "Mono", "run_id", "gonl_id")])

metadata$uuid <- colnames(betas)
head(metadata)
metadata <- merge(metadata, idat[,c("idat", "run_id")], by="run_id")
str(metadata)

metadata$set <- "unrelated"

write.table(metadata, "EGA_upload_DNAm.csv", row.names=FALSE, quote=FALSE, sep="\t")

##EGA_upload_unrelated_RNA
data(rnaSeqData_ReadCounts_BIOS_cleaned)

##for fastq path and  insert size
fastq <- getView("getFastq", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB)
fastq <- fastq[!duplicated(fastq$run_id),]
dim(fastq)

fastq$fastq <- apply(fastq[, c("R1", "R2")], 1, paste, collapse=", ")

bam <- getView("getBAM", usrpwd=RP3_MDB_USRPWD, url=RP3_RDB)
bam <- bam[grepl("v2.1.3", bam$path),] ##freeze 2 alignments

run <- merge(fastq[, c("run_id", "fastq", "insert_size")], bam[, c("run_id", "path")], by="run_id")
colnames(run)[4] <- "bam"
head(run)

dim(run)

metadata <- as.data.frame(colData(counts)[, c("Sex", "Sampling_Age", "Baso_Perc", "Eos_Perc", "Neut_Perc", "LUC_Perc", "Lymph_Perc", "Mono_Perc", "Baso", "Eos", "Neut", "LUC", "Lymph", "Mono", "run_id", "gonl_id")])

metadata$uuid <- colnames(counts)
metadata[1:5,]
metadata <- merge(metadata, run, by="run_id")
metadata[1:5,]
dim(metadata)

##some how this run escaped should not be merged
setdiff(metadata$run_id, run$run_id)
subset(metadata, run_id == "BD1P0GACXX-4-4_BC52YAACXX-8-4")
ids <- getView("getIds", usrpwd=RP3_MDB_USRPWD, url=RP3_MDB)
subset(ids, uuid == "BIOSB6925833")

metadata$set <- "unrelated"
colnames(metadata)

##EGA_upload_GoNL_RNA
data(rnaSeqData_ReadCounts_GoNL)

md <- as.data.frame(colData(counts)[, c("Sex", "Sampling_Age", "Baso_Perc", "Eos_Perc", "Neut_Perc", "LUC_Perc", "Lymph_Perc", "Mono_Perc", "Baso", "Eos", "Neut", "LUC", "Lymph", "Mono", "run_id", "gonl_id")])

md$uuid <- colnames(counts)

md <- merge(md, run, by="run_id")
md[1:5,]

md$set <- "GoNL"

colnames(md)

metadata <- rbind(metadata, md)

dim(metadata)

write.table(metadata, "EGA_upload_RNA.csv", row.names=FALSE, quote=FALSE, sep="\t")

