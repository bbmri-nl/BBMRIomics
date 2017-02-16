#!/bin/bash

# This script downloads impute2 files from the SRM
# and converts them to tabix indexed files
# 2015-02-10, M. van Iterson <m.van_iterson@lumc.nl>

URL="srm://srm.grid.sara.nl/pnfs/grid.sara.nl/data/bbmri.nl/RP3/GWAS_ImputationGoNLv5"
BIOBANK=( PAN LL LLS RS CODAM NTR_B37_Aff6_Epigen NTR_B37_AC_Epigen_CAUT )
DATE=( 20140306 20140306 20140402 20140306 20140306 20140606 20140606 )
MAF=0.05
INFO=0.4

( 

##install tabix
##git clone https://github.com/samtools/tabix.git
##cd tabix
##make
##cd ..

for i in `seq 1 2 3 4 5 6 7`;
do
    echo "Biobank ${BIOBANK[i]}"
    if [[ "${BIOBANK[i]}" =~ "NTR" ]] ; then
         INSFILE="${URL}/NTR/${BIOBANK[i]}-imputed-${DATE[i]}.sample"
     else
         INSFILE="${URL}/${BIOBANK[i]}/${BIOBANK[i]}-imputed-${DATE[i]}.sample"
    fi
    OUTSFILE="$( basename $INSFILE)"
    srmcp --server_mode=passive $INSFILE "file:///$OUTSFILE"

    for CHR in `seq 1 22`;
    do
        echo "Chromosome ${CHR}"
        if [[ "${BIOBANK[i]}" =~ "NTR" ]] ; then
            INGFILE="${URL}/NTR/${BIOBANK[i]}-imputed-chr$CHR-${DATE[i]}"
        else
            INGFILE="${URL}/${BIOBANK[i]}/${BIOBANK[i]}-imputed-chr$CHR-${DATE[i]}"
        fi
	OUTGFILE="$( basename $INGFILE)"

        echo "copy file: ${INGFILE}"      
        echo "to file: ${OUTGFILE}"      
	
	srmcp --server_mode=passive $INGFILE "file:///$OUTGFILE"                    #using srmcp from lsg UI
	srmcp --server_mode=passive ${INGFILE}".md5" "file:///${OUTGFILE}.md5"
	srmcp --server_mode=passive ${INGFILE}"_info" "file:///${OUTGFILE}_info"
	srmcp --server_mode=passive ${INGFILE}"_info.md5" "file:///${OUTGFILE}_info.md5"

	MD5ORIG=$( cut -d ' ' -f 1 ${OUTGFILE}".md5")
	MD5COPY=$(md5sum $OUTGFILE | cut -d ' ' -f 1)	
      
	if [ "$MD5ORIG" != "$MD5COPY" ]
	then
	    echo "MD5 mismatch on file: $OUTGFILE"  ##redo copying?
	fi

	MD5ORIG=$( cut -d ' ' -f 1 ${OUTGFILE}"_info.md5")
	MD5COPY=$(md5sum ${OUTGFILE}"_info" | cut -d ' ' -f 1)	
      
	if [ "$MD5ORIG" != "$MD5COPY" ]
	then
	    echo "MD5 mismatch on file: ${OUTGFILE}_info"  ##redo copying?
	fi

	awk -v chr="$CHR" '{
            printf chr"\t"$2"\t"$3"\t"$4"\t"$5; 
            for(i=6; i<NF; i+=3) {
                 printf "\t"$(i+0)*0+$(i+1)*1+$(i+2)*2
             }; 
            printf "\n"}' $OUTGFILE > tmp                   
	
	##add sample info
        tail -n +3 $OUTSFILE | cut -d ' ' -f2 > header
	echo -e 'chr\nrsid\npos\nref\nalt' > header1
	cat header >> header1
	tr '\n' '\t' < header1 > header2
	sed -i -e '$a\' header2
	sed -i 's/\t$//g' header2

	cat header2 tmp > tmp1
	tr ' ' '\t' < ${OUTGFILE}"_info" > ${OUTGFILE}"_info.tmp"

        paste ${OUTGFILE}"_info.tmp" tmp1 > $OUTGFILE ##combine info and genotypes
	
	awk -v info="$INFO" -F' ' '{if($5>info) print}' $OUTGFILE > tmp ##filter on info	

        sort -n -k 10 tmp > $OUTGFILE                                               
        /home/miterson/tabix/bgzip $OUTGFILE
        /home/miterson/tabix/tabix ${OUTGFILE}".gz" -s 8 -b 10 -e 10 -S 1
	
	rm *.md5 tmp* header* *info* 
    done
done

) > impute2tabix.log 2>&1
