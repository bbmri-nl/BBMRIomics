#Makefile for linking datasets and extract data from databases
#CONFIG_FILE should contain update urls and paths

ROOT_DIR:=$(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
CONFIG_FILE = ${ROOT_DIR}/bbmriomics.conf

BIOBANKS = RS PAN CODAM LLS LL NTR

VM_BASE_DATA = $(shell grep "^VM_BASE_DATA: " ${CONFIG_FILE} | cut -d" " -f2 | sed -e 's/^[ \t]*//')
RP3_MDB = $(shell grep "^RP3_MDB: " ${CONFIG_FILE} | cut -d" " -f2 | sed -e 's/^[ \t]*//')
RP3_RDB = $(shell grep "^RP3_RDB: " ${CONFIG_FILE} | cut -d" " -f2 | sed -e 's/^[ \t]*//')

R_PACKAGE_DIR = ${R_LIBS_SITE}/BBMRIomics

#all: init rp3 rp4 views
all: init rp3 rp4

init:
	mkdir -p $(R_PACKAGE_DIR)/views
	echo "VM_BASE_DATA=" $(VM_BASE_DATA)
	echo "RP3_MDB=" $(RP3_MDB)
	echo "RP3_RDB=" $(RP3_RDB)

rp3: $(BIOBANKS) BIOS

$(BIOBANKS): 
	ln -fs $(VM_BASE_DATA)/IlluminaHumanMethylation450k/450k/$@/methData_Betas_$@_F2_cleaned.RData $(R_PACKAGE_DIR)/data/methData_Betas_$@_F2_cleaned.RData
	ln -fs $(VM_BASE_DATA)/IlluminaHumanMethylation450k/450k/$@/methData_Mvalues_$@_F2_cleaned.RData $(R_PACKAGE_DIR)/data/methData_Mvalues_$@_F2_cleaned.RData
	ln -fs $(VM_BASE_DATA)/RNASeq/v2.1.3/gene_read/rnaSeqData_ReadCounts_$@_cleaned.RData $(R_PACKAGE_DIR)/data/rnaSeqData_ReadCounts_$@_cleaned.RData
	ln -fs $(VM_BASE_DATA)/IlluminaHumanMethylation450k/450k/$@/methData_Betas_$@_Freeze2_unrelated.RData $(R_PACKAGE_DIR)/data/methData_Betas_$@_Freeze2_unrelated.RData
	ln -fs $(VM_BASE_DATA)/IlluminaHumanMethylation450k/450k/$@/methData_Mvalues_$@_Freeze2_unrelated.RData $(R_PACKAGE_DIR)/data/methData_Mvalues_$@_Freeze2_unrelated.RData
	ln -fs $(VM_BASE_DATA)/RNASeq/v2.1.3/gene_read/rnaSeqData_ReadCounts_$@_Freeze2_unrelated.RData $(R_PACKAGE_DIR)/data/rnaSeqData_ReadCounts_$@_Freeze2_unrelated.RData
	ln -fs $(VM_BASE_DATA)/RNASeq/GRCh38/gene_read/rnaSeqData_ReadCounts_$@_Freeze2_unrelated_GRCh38.RData $(R_PACKAGE_DIR)/data/rnaSeqData_ReadCounts_$@_Freeze2_unrelated_GRCh38.RData

BIOS:
	ln -fs $(VM_BASE_DATA)/IlluminaHumanMethylation450k/450k/methData_Betas_BIOS_F2_cleaned.RData $(R_PACKAGE_DIR)/data/methData_Betas_BIOS_F2_cleaned.RData
	ln -fs $(VM_BASE_DATA)/IlluminaHumanMethylation450k/450k/methData_Mvalues_BIOS_F2_cleaned.RData $(R_PACKAGE_DIR)/data/methData_Mvalues_BIOS_F2_cleaned.RData
	ln -fs $(VM_BASE_DATA)/RNASeq/v2.1.3/gene_read/rnaSeqData_ReadCounts_BIOS_cleaned.RData $(R_PACKAGE_DIR)/data/rnaSeqData_ReadCounts_BIOS_cleaned.RData
	ln -fs $(VM_BASE_DATA)/RNASeq/v2.1.3/gene_read/rnaSeqData_ReadCounts_GoNL.RData $(R_PACKAGE_DIR)/data/rnaSeqData_ReadCounts_GoNL.RData
	ln -fs $(VM_BASE_DATA)/IlluminaHumanMethylation450k/450k/methData_Betas_BIOS_Freeze2_unrelated.RData $(R_PACKAGE_DIR)/data/methData_Betas_BIOS_Freeze2_unrelated.RData
	ln -fs $(VM_BASE_DATA)/IlluminaHumanMethylation450k/450k/methData_Mvalues_BIOS_Freeze2_unrelated.RData $(R_PACKAGE_DIR)/data/methData_Mvalues_BIOS_Freeze2_unrelated.RData
	ln -fs $(VM_BASE_DATA)/RNASeq/v2.1.3/gene_read/rnaSeqData_ReadCounts_BIOS_Freeze2_unrelated.RData $(R_PACKAGE_DIR)/data/rnaSeqData_ReadCounts_BIOS_Freeze2_unrelated.RData
	ln -fs $(VM_BASE_DATA)/RNASeq/v2.1.3/gene_read/rnaSeqData_ReadCounts_BIOS_Freeze2_GoNL.RData $(R_PACKAGE_DIR)/data/rnaSeqData_ReadCounts_BIOS_Freeze2_GoNL.RData
	ln -fs $(VM_BASE_DATA)/RNASeq/GRCh38/gene_read/rnaSeqData_ReadCounts_BIOS_Freeze2_unrelated_GRCh38.RData $(R_PACKAGE_DIR)/data/rnaSeqData_ReadCounts_BIOS_Freeze2_unrelated_GRCh38.RData
	ln -fs $(VM_BASE_DATA)/RNASeq/GRCh38/gene_read/rnaSeqData_ReadCounts_BIOS_Freeze2_GoNL_GRCh38.RData $(R_PACKAGE_DIR)/data/rnaSeqData_ReadCounts_BIOS_Freeze2_GoNL_GRCh38.RData

rp4:
	ln -fs $(VM_BASE_DATA)/RP4/metabolomics_RP3RP4_overlap.RData $(R_PACKAGE_DIR)/data/metabolomics_RP3RP4_overlap.RData

#views:
#	curl -X GET $(RP3_MDB)/_design/EGA/_view/freeze1RNASeq?reduce=false -u $(USRPWD)  -k -g >  $(R_PACKAGE_DIR)/views/EGA_freeze1RNASeq.json
#	curl -X GET $(RP3_MDB)/_design/EGA/_view/freeze1Methylation?reduce=false -u $(USRPWD)  -k -g >  $(R_PACKAGE_DIR)/views/EGA_freeze1Methylation.json
#	curl -X GET $(RP3_MDB)/_design/EGA/_view/freeze2RNASeq?reduce=false -u $(USRPWD)  -k -g >  $(R_PACKAGE_DIR)/views/EGA_freeze2RNASeq.json
#	curl -X GET $(RP3_MDB)/_design/EGA/_view/freeze2Methylation?reduce=false -u $(USRPWD)  -k -g >  $(R_PACKAGE_DIR)/views/EGA_freeze2Methylation.json
#	curl -X GET $(RP3_MDB)/_design/Files/_view/getFastq?reduce=false -u $(USRPWD)  -k -g >  $(R_PACKAGE_DIR)/views/Files_getFastq.json
#	curl -X GET $(RP3_MDB)/_design/Files/_view/getIdat?reduce=false -u $(USRPWD)  -k -g >  $(R_PACKAGE_DIR)/views/Files_getIdat.json
#	curl -X GET $(RP3_MDB)/_design/Identifiers/_view/getIds?reduce=false -u $(USRPWD)  -k -g >  $(R_PACKAGE_DIR)/views/Identifiers_getIds.json
#	curl -X GET $(RP3_MDB)/_design/Phenotypes/_view/allPhenotypes?reduce=false -u $(USRPWD)  -k -g >  $(R_PACKAGE_DIR)/views/Phenotypes_allPhenotypes.json
#	curl -X GET $(RP3_MDB)/_design/Phenotypes/_view/cellCounts?reduce=false -u $(USRPWD)  -k -g >  $(R_PACKAGE_DIR)/views/Phenotypes_cellCounts.json
#	curl -X GET $(RP3_MDB)/_design/Phenotypes/_view/minimalPhenotypes?reduce=false -u $(USRPWD)  -k -g >  $(R_PACKAGE_DIR)/views/Phenotypes_minimalPhenotypes.json
#	curl -X GET $(RP3_MDB)/_design/Runs/_view/getGenotypes?reduce=false -u $(USRPWD)  -k -g >  $(R_PACKAGE_DIR)/views/Runs_getGenotypes.json
#	curl -X GET $(RP3_MDB)/_design/Runs/_view/getMethylationRuns?reduce=false -u $(USRPWD)  -k -g >  $(R_PACKAGE_DIR)/views/Runs_getMethylationRuns.json
#	curl -X GET $(RP3_MDB)/_design/Runs/_view/getRNASeqRuns?reduce=false -u $(USRPWD)  -k -g >  $(R_PACKAGE_DIR)/views/Runs_getRNASeqRuns.json
#	curl -X GET $(RP3_MDB)/_design/Runs/_view/getImputations?reduce=false -u $(USRPWD)  -k -g >  $(R_PACKAGE_DIR)/views/Runs_getImputations.json
#	curl -X GET $(RP3_MDB)/_design/Runs/_view/getRelations?reduce=false -u $(USRPWD)  -k -g >  $(R_PACKAGE_DIR)/views/Runs_getRelations.json
#	curl -X GET $(RP3_MDB)/_design/Samplesheets/_view/rnaseqSamplesheet?reduce=false -u $(USRPWD)  -k -g >  $(R_PACKAGE_DIR)/views/Samplesheets_rnaseqSamplesheet.json
#	curl -X GET $(RP3_MDB)/_design/Samplesheets/_view/methylationSamplesheet?reduce=false -u $(USRPWD)  -k -g >  $(R_PACKAGE_DIR)/views/Samplesheets_methylationSamplesheet.json
#	curl -X GET $(RP3_MDB)/_design/Verification/_view/md5sum?reduce=false -u $(USRPWD)  -k -g >  $(R_PACKAGE_DIR)/views/Verification_md5sum.json
#	curl -X GET $(RP3_RDB)/_design/stats/_view/getStats?reduce=false -u $(USRPWD)  -k -g >  $(R_PACKAGE_DIR)/views/stats_getStats.json
#	curl -X GET $(RP3_RDB)/_design/stats/_view/getBAM?reduce=false -u $(USRPWD)  -k -g >  $(R_PACKAGE_DIR)/views/stats_getBAM.json

stop:
	@echo ** [BBMRIomics] Installing without linking data and generating stored-views

