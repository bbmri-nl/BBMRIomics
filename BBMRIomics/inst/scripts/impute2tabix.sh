#!/bin/bash

# This script downloads impute2 files from the SRM
# and converts them to tabix indexed files
# 2015-02-10, M. van Iterson <m.van_iterson@lumc.nl>

URL="srm://srm.grid.sara.nl/pnfs/grid.sara.nl/data/bbmri.nl/RP3/GWAS_ImputationGoNLv5"
BIOBANK=( PAN LL LLS RS CODAM NTR_B37_Aff6_Epigen NTR_B37_AC_Epigen_CAUT )
DATE=( 20140306 20140306 20140402 20140306 20140306 20140606 20140606 )
MAF=0.05
CER=0.4

#Certificate. Needed to download bam files from srm storage.  
export X509_USER_PROXY=/home/miterson/x509up_u34722

( 

##install tabix
##git clone https://github.com/samtools/tabix.git
##cd tabix
##make
##cd ..

for i in `seq 1 7`;
do
    echo "Biobank ${BIOBANK[i]}"
    for CHR in `seq 1 22`;
    do
        echo "Chromosome ${CHR}"
        if [[ "${BIOBANK[i]}" =~ "NTR" ]] ; then
            INFILE="${URL}/NTR/${BIOBANK[i]}-imputed-chr$CHR-${DATE[i]}"
        else
            INFILE="${URL}/${BIOBANK[i]}/${BIOBANK[i]}-imputed-chr$CHR-${DATE[i]}"
        fi
	OUTFILE="$( basename $INFILE)"

        echo "copy file: ${INFILE}"      
        echo "to file: ${OUTFILE}"      
	
	srmcp --server_mode=passive $INFILE "file:///$OUTFILE"                    #using srmcp from lsg UI
	srmcp --server_mode=passive ${INFILE}".md5" "file:///${OUTFILE}.md5"
	srmcp --server_mode=passive ${INFILE}"_info" "file:///${OUTFILE}_info"
	srmcp --server_mode=passive ${INFILE}"_info.md5" "file:///${OUTFILE}_info.md5"

	MD5ORIG=$( cut -d ' ' -f 1 ${OUTFILE}".md5")
	MD5COPY=$(md5sum $OUTFILE | cut -d ' ' -f 1)	
      
	if [ "$MD5ORIG" != "$MD5COPY" ]
	then
	    echo "MD5 mismatch on file: $FILE"  ##redo copying?
	fi

	MD5ORIG=$( cut -d ' ' -f 1 ${OUTFILE}"_info.md5")
	MD5COPY=$(md5sum ${OUTFILE}"_info" | cut -d ' ' -f 1)	
      
	if [ "$MD5ORIG" != "$MD5COPY" ]
	then
	    echo "MD5 mismatch on file: ${OUTFILE}_info"  ##redo copying?
	fi

        paste ${OUTFILE}"_info" $OUTFILE > tmp ##combine info and genotypes
	
	##awk -v maf="$MAF" -v cer="$CER" -F' ' '{if($4>maf && $4 < 1-maf && $6>cer) print}' tmp > ${OUTFILE}	

	awk -v chr="$CHR" '{
            printf chr"\t"$2"\t"$3"\t"$4"\t"$5; 
            for(i=6; i<NF; i+=3) {
                 printf "\t"$(i+0)*0+$(i+1)*1+$(i+2)*2
             }; 
            printf "\n"}' tmp > $OUTFILE                   


        ##sort -n -k 3 tmp > $OUTFILE                                               
        /home/miterson/tabix/bgzip $OUTFILE
        /home/miterson/tabix/tabix ${OUTFILE}".gz" -s 1 -b 3 -e 3
	
	rm *.md5 tmp
    done
done

) > impute2tabix.log 2>&1
