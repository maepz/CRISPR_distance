#!/bin/bash
path=$1
cd ${path}

declare -a arr=("PleD" "LpxA" "TufB")
for locus in "${arr[@]}"
do


i=1
for file in `ls ${locus}/*merged.fastq.gz`
do

FILE="${file/%_merged.fastq.gz/}"
SAMPLE=`echo $FILE | cut -d'/' -f2`
LIB=`echo $SAMPLE | cut -d'_' -f3`
REGION=`echo $LIB | cut -d'-' -f1 | cut -d'_' -f2`
IND=`echo $LIB | cut -d'_' -f2`
#echo $FILE
echo $SAMPLE
#echo $LIB
#echo $REGION
#echo $IND
#echo $i
i=$((i + 1))

# remap best reads
bowtie2 -x ../bowtie2_indices/Symb_1 -U ${FILE}_merged.fastq.gz --very-sensitive-local -qseq -X 600 --rg-id ${i} --rg LB:${LIB} --rg SM:${LIB} --rg PI:275 --rg CN:${IND} | samtools view -bS - > ${FILE}.bam


samtools sort ${FILE}.bam -o ${FILE}.sorted.bam
samtools index ${FILE}.sorted.bam
rm ${FILE}.bam

done


ls ${locus}/${locus}*.sorted.bam > ${locus}/bamfileslist.txt
samtools merge -f -b ${locus}/bamfileslist.txt ${locus}/ALL_${locus}.sorted.bam
samtools index ${locus}/ALL_${locus}.sorted.bam

samtools mpileup -I -f ../bowtie2_indices/Symb_1.fasta ${locus}/ALL_${locus}.sorted.bam | java -jar ../VarScan.v2.3.9.jar mpileup2snp --min-coverage 100 --min-reads2 10 --min-avg-qual 25 --min-var-freq 0.01 --output-vcf 1 > ALL_${locus}.vcf

done