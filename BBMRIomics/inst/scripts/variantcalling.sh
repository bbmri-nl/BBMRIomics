IFS=$'\n' read -d '' -r -a SAMPLES < $1

PROJDIR=/exports/molepi/users/mvaniterson/variantcalling
URL=https://fly1.grid.sara.nl:2882/pnfs/grid.sara.nl/data/bbmri.nl/RP3/RNASeq//v2.1.3/

MYPROXY="${PROJDIR}/x509up_u34722"
CURL="curl --CApath /etc/grid-security/certificates/ -E ${MYPROXY} -L -k"
CONFIG="${PROJDIR}/config.yml"
PICARD="java -Xmx4G -jar /usr/local/sasc/programs/picard-tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups"

REFFASTA=/exports/molepi/RNA_seq_Leen/project-152-LeenRNAseq/data/reference/reference.fasta
BEDFILE=/exports/molepi/users/mvaniterson/variantcalling/merge.sorted.merge.bed

##update proxy first run startGridsession from ui and copy to VM
# scp bios16:/tmp/x509up_u34722 .

##run the variant calling
##download bams and do the variantcalling
mkdir -p ${PROJDIR}/output

for SAMPLE in "${SAMPLES[@]}"
do
    SCRIPT=${PROJDIR}/${SAMPLE}.sh
    BAM="${SAMPLE}.mdup.bam"
    BAMIN="${URL}/${BAM}"
    
cat > ${SCRIPT} << EOT
#!/bin/bash
#$ -S /bin/bash
#$ -q all.q
#$ -N VariantCalling
#$ -l h_vmem=8G
#$ -wd ${PROJDIR}/output
#$ -j Y
#$ -V
#$ -m be
#$ -M m.van_iterson@lumc.nl


echo ${SAMPLE} > ${SAMPLE}.txt

echo "...downloading"
${CURL} ${BAMIN} | \
samtools mpileup -f ${REFFASTA} -l ${BEDFILE} -d 1000000 -s -B - | java -XX:+UseParallelOldGC -XX:ParallelGCThreads=4 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx4G -jar /usr/local/VarScan/VarScan.v2.3.7.jar mpileup2cns -strand-filter 0 --output-vcf 1 -vcf-sample-list ${SAMPLE}.txt | bgzip -c > ${SAMPLE}.vcf.gz
tabix -p vcf ${SAMPLE}.vcf.gz
rm ${SAMPLE}.txt
EOT
    
    qsub ${SCRIPT}
    
done

