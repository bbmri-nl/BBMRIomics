#!/bin/bash

# Genotyping
#
# M. van Iterson <m.van_iterson@lumc.nl>
RSCRIPT=/usr/lib64/microsoft-r/3.3/lib64/R/bin/Rscript
OUT=/virdir/Backup/RP3_analysis/SwapDetection/

##intra
$RSCRIPT genotyping.R --typex HRC --cohort ALL --out ${OUT} --filex ${OUT}/HighQualPositions.GCRh37.bed
$RSCRIPT genotyping.R --typex GoNL --cohort ALL --out ${OUT} --filex ${OUT}/HighQualPositions.GCRh37.bed
$RSCRIPT genotyping.R --typex DNAm --cohort ALL --out ${OUT}
$RSCRIPT genotyping.R --typex RNA --cohort ALL --out ${OUT} --filex ${OUT}/output_RNA_2600.vcf

##inter
$RSCRIPT genotyping.R --typex RNA --typey HRC --cohort ALL --out ${OUT} --filex ${OUT}/output_RNA_2600.vcf
$RSCRIPT genotyping.R --typex RNA --typey GoNL --cohort ALL --out ${OUT} --filex ${OUT}/output_RNA_2600.vcf
$RSCRIPT genotyping.R --typex DNAm --typey HRC --cohort ALL --out ${OUT}
$RSCRIPT genotyping.R --typex DNAm --typey GoNL --cohort ALL --out ${OUT}
$RSCRIPT genotyping.R --typex HRC --typey GoNL --cohort ALL --filex ${OUT}/HighQualPositions.GCRh37.bed --out ${OUT}
