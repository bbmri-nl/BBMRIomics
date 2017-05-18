
SRC=/virdir/Scratch/RP3_data
BIOBANKS=RS PAN CODAM LLS LL NTR

USRPWD=''

MDB=https://metadatabase.bbmrirp3-lumc.surf-hosted.nl/bios
RDB=https://metadatabase.bbmrirp3-lumc.surf-hosted.nl/rp3_analysis

R_PACKAGE_DIR=/opt/R/microsoft-r/3.3/lib64/R/library/BBMRIomics

all: init rp3 rp4 views

init:
	mkdir -p $(R_PACKAGE_DIR)/views

rp3: $(BIOBANKS) BIOS

$(BIOBANKS): 
	ln -fs $(SRC)/IlluminaHumanMethylation450k/450k/$@/methData_Betas_$@_F2_cleaned.RData $(R_PACKAGE_DIR)/data/methData_Betas_$@_F2_cleaned.RData
	ln -fs $(SRC)/IlluminaHumanMethylation450k/450k/$@/methData_Mvalues_$@_F2_cleaned.RData $(R_PACKAGE_DIR)/data/methData_Mvalues_$@_F2_cleaned.RData
	ln -fs $(SRC)/RNASeq/v2.1.3/gene_read/rnaSeqData_ReadCounts_$@_cleaned.RData $(R_PACKAGE_DIR)/data/rnaSeqData_ReadCounts_$@_cleaned.RData

BIOS:
	ln -fs $(SRC)/IlluminaHumanMethylation450k/450k/methData_Betas_BIOS_F2_cleaned.RData $(R_PACKAGE_DIR)/data/methData_Betas_BIOS_F2_cleaned.RData
	ln -fs $(SRC)/IlluminaHumanMethylation450k/450k/methData_Mvalues_BIOS_F2_cleaned.RData $(R_PACKAGE_DIR)/data/methData_Mvalues_BIOS_F2_cleaned.RData
	ln -fs $(SRC)/RNASeq/v2.1.3/gene_read/rnaSeqData_ReadCounts_BIOS_cleaned.RData $(R_PACKAGE_DIR)/data/rnaSeqData_ReadCounts_BIOS_cleaned.RData

rp4:
	ln -fs $(SRC)/RP4/metabolomics_RP3RP4_overlap.RData $(R_PACKAGE_DIR)/data/metabolomics_RP3RP4_overlap.RData

views:
	curl -X GET $(MDB)/_design/EGA/_view/freeze1RNASeq?reduce=false -u $(USRPWD)  -k -g >  $(R_PACKAGE_DIR)/views/EGA_freeze1RNASeq.json
	curl -X GET $(MDB)/_design/EGA/_view/freeze1Methylation?reduce=false -u $(USRPWD)  -k -g >  $(R_PACKAGE_DIR)/views/EGA_freeze1Methylation.json
	curl -X GET $(MDB)/_design/EGA/_view/freeze2RNASeq?reduce=false -u $(USRPWD)  -k -g >  $(R_PACKAGE_DIR)/views/EGA_freeze2RNASeq.json
	curl -X GET $(MDB)/_design/EGA/_view/freeze2Methylation?reduce=false -u $(USRPWD)  -k -g >  $(R_PACKAGE_DIR)/views/EGA_freeze2Methylation.json
	curl -X GET $(MDB)/_design/Files/_view/getFastq?reduce=false -u $(USRPWD)  -k -g >  $(R_PACKAGE_DIR)/views/Files_getFastq.json
	curl -X GET $(MDB)/_design/Files/_view/getIdat?reduce=false -u $(USRPWD)  -k -g >  $(R_PACKAGE_DIR)/views/Files_getIdat.json
	curl -X GET $(MDB)/_design/Identifiers/_view/getIds?reduce=false -u $(USRPWD)  -k -g >  $(R_PACKAGE_DIR)/views/Identifiers_getIds.json
	curl -X GET $(MDB)/_design/Phenotypes/_view/allPhenotypes?reduce=false -u $(USRPWD)  -k -g >  $(R_PACKAGE_DIR)/views/Phenotypes_allPhenotypes.json
	curl -X GET $(MDB)/_design/Phenotypes/_view/cellCounts?reduce=false -u $(USRPWD)  -k -g >  $(R_PACKAGE_DIR)/views/Phenotypes_cellCounts.json
	curl -X GET $(MDB)/_design/Phenotypes/_view/minimalPhenotypes?reduce=false -u $(USRPWD)  -k -g >  $(R_PACKAGE_DIR)/views/Phenotypes_minimalPhenotypes.json
	curl -X GET $(MDB)/_design/Runs/_view/getGenotypes?reduce=false -u $(USRPWD)  -k -g >  $(R_PACKAGE_DIR)/views/Runs_getGenotypes.json
	curl -X GET $(MDB)/_design/Runs/_view/getMethylationRuns?reduce=false -u $(USRPWD)  -k -g >  $(R_PACKAGE_DIR)/views/Runs_getMethylationRuns.json
	curl -X GET $(MDB)/_design/Runs/_view/getRNASeqRuns?reduce=false -u $(USRPWD)  -k -g >  $(R_PACKAGE_DIR)/views/Runs_getRNASeqRuns.json
	curl -X GET $(MDB)/_design/Runs/_view/getImputations?reduce=false -u $(USRPWD)  -k -g >  $(R_PACKAGE_DIR)/views/Runs_getImputations.json
	curl -X GET $(MDB)/_design/Runs/_view/getRelations?reduce=false -u $(USRPWD)  -k -g >  $(R_PACKAGE_DIR)/views/Runs_getRelations.json
	curl -X GET $(MDB)/_design/Samplesheets/_view/rnaseqSamplesheet?reduce=false -u $(USRPWD)  -k -g >  $(R_PACKAGE_DIR)/views/Samplesheets_rnaseqSamplesheet.json
	curl -X GET $(MDB)/_design/Samplesheets/_view/methylationSamplesheet?reduce=false -u $(USRPWD)  -k -g >  $(R_PACKAGE_DIR)/views/Samplesheets_methylationSamplesheet.json
	curl -X GET $(MDB)/_design/Verification/_view/md5sum?reduce=false -u $(USRPWD)  -k -g >  $(R_PACKAGE_DIR)/views/Verification_md5sum.json
	curl -X GET $(RDB)/_design/stats/_view/getStats?reduce=false -u $(USRPWD)  -k -g >  $(R_PACKAGE_DIR)/views/stats_getStats.json
	curl -X GET $(RDB)/_design/stats/_view/getBAM?reduce=false -u $(USRPWD)  -k -g >  $(R_PACKAGE_DIR)/views/stats_getBAM.json

stop:
	@echo ** [BBMRIomics] Installing without linking data and generating stored-views

