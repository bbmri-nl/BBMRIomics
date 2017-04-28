#!/bin/bash

# Genotyping
# 
# 2017-26-01, M. van Iterson <m.van_iterson@lumc.nl>
RSCRIPT=Rscript
COHORTS=(PAN LL LLS RS CODAM NTR)
#TYPES=(DNAm RNA HRC GoNL)
TYPES=(DNAm)

echo $RSCRIPT --version

#inter
# for cohort in ${COHORTS}
# do
#     $RSCRIPT genotyping.R --typex DNAm --typey HRC --cohort $cohort
#     $RSCRIPT genotyping.R --typex DNAm --typey GoNL --cohort $cohort
#     $RSCRIPT genotyping.R --typex RNA --typey HRC --cohort $cohort
#     $RSCRIPT genotyping.R --typex RNA --typey GoNL --cohort $cohort
# done

#intra
for type in ${TYPES}
do
    $RSCRIPT genotyping.R --typex $type
done





