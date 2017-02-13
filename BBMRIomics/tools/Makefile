
SRC=/virdir/Scratch/RP3_data
BIOBANKS=RS PAN CODAM LLS LL NTR
MDB=https://metadatabase.bbmrirp3-lumc.surf-hosted.nl:6984/bios
RDB=https://metadatabase.bbmrirp3-lumc.surf-hosted.nl:6984/rp3_analysis
USRPWD='usrpwd'
R_PACKAGE_DIR=/usr/lib64/microsoft-r/3.3/lib64/R/library/BBMRIomics

ifeq ($(shell hostname),packer-ubuntu-16)
	run=init methylation
##rnaseq metabolomics views
else
	run=stop
endif

all: $(run)

init:
	mkdir -p $(R_PACKAGE_DIR)/views

methylation: $(BIOBANKS) BIOS

$(BIOBANKS): 
	ln -s $(SRC)/IlluminaHumanMethylation450k/450k/$@/methData_Betas_$@_F2.RData $(R_PACKAGE_DIR)/data/methData_Betas_$@_F2_cleaned.RData
	ln -s $(SRC)/IlluminaHumanMethylation450k/450k/$@/methData_Mvalues_$@_F2.RData $(R_PACKAGE_DIR)/data/methData_Mvalues_$@_F2_cleaned.RData

BIOS:
#	ln -s $(SRC)/IlluminaHumanMethylation450k/450k/methData_Betas_BIOS_F2.RData $(R_PACKAGE_DIR)/data/methData_Betas_BIOS_F2.RData
#	ln -s $(SRC)/IlluminaHumanMethylation450k/450k/methData_Mvalues_BIOS_F2.RData $(R_PACKAGE_DIR)/data/methData_Mvalues_BIOS_F2.RData

rnaseq: genes exons

genes:
	ln -s $(SRC)/RNASeq/run_01/genecounts/rnaSeqData_freeze1_06032015BIOS.RData $(R_PACKAGE_DIR)/data/rnaSeqData_freeze1_06032015BIOS.RData
	ln -s $(SRC)/RNASeq/v2.1.3/gene_base/geneCountsF2.RData $(R_PACKAGE_DIR)/data/geneCountsF2.RData

exons:
	ln -s $(SRC)/RNASeq/run_01/exoncounts/rnaSeqData_freeze1_exon_14042015BIOS.RData $(R_PACKAGE_DIR)/data/rnaSeqData_freeze1_exon_14042015BIOS.RData
	ln -s $(SRC)/RNASeq/v2.1.3/exon_base/exonCountsF2.RData $(R_PACKAGE_DIR)/data/exonCountsF2.RData

metabolomics:
	ln -s $(SRC)/RP4/metabolomics_RP3RP4_overlap.RData $(R_PACKAGE_DIR)/data/metabolomics_RP3RP4_overlap.RData

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

stop:
	@echo ** [BBMRIomics] Installing without linking data and generating stored-views

