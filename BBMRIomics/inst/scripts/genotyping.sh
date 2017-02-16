#!/bin/bash

# Genotyping
# 
# 2017-26-01, M. van Iterson <m.van_iterson@lumc.nl>
RSCRIPT=/usr/local/R/R-3.2.0/bin/Rscript
COHORTS=(PAN LL LLS RS CODAM NTR)
TYPES=(DNAm RNA HRC GoNL)

#inter
for cohort in ${COHORTS}
do
    $RSCRIPT genotyping.R -x DNAm -y HRC -c $cohort
    $RSCRIPT genotyping.R -x DNAm -y GoNL -c $cohort
    $RSCRIPT genotyping.R -x RNA -y HRC -c $cohort
    $RSCRIPT genotyping.R -x RNA -y GoNL -c $cohort
done

#intra
for type in ${TYPES}
do
    $RSCRIPT genotyping.R -x $type
done




