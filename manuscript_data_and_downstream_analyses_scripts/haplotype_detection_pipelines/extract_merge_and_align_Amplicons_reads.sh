#!/bin/bash
DIR=Amplicons
DIR=$1

cd $DIR

mkdir V4
mkdir LpxA
mkdir PleD
mkdir TufB
mkdir unmapped


for file in `ls *.bam`
do
SAMPLE="${file/%.bam/}"
echo " Working of sample ${SAMPLE} ..."


#############
echo "
extracting unmapped reads ...
The dependencies required here are: samtools (http://www.htslib.org/download/), and bedtools (https://bedtools.readthedocs.io/en/latest/index.html)"
#############

samtools view -u -f 12 -F 256  ${SAMPLE}.bam > unmapped/unmapped.${SAMPLE}.bam
bamToFastq -i unmapped/unmapped.${SAMPLE}.bam -fq unmapped/unmapped.${SAMPLE}_R1.fastq -fq2 unmapped/unmapped.${SAMPLE}_R2.fastq
gzip unmapped/unmapped.${SAMPLE}_R1.fastq
gzip unmapped/unmapped.${SAMPLE}_R2.fastq
rm unmapped/unmapped.${SAMPLE}.bam

#############
echo "
extracting reads that mapped to reference genome. The extracted reads are only those where both paired-end mapped "
#############

samtools index ${SAMPLE}.bam


##########
echo "extracting reads that mapped V4 ...."
###########


samtools view -h -f 2 -F 4 ${SAMPLE}.bam Ga0074115_140:2462-3003 > V4/V4_${SAMPLE}.sorted.bam
samtools sort -n V4/V4_${SAMPLE}.sorted.bam -o V4/V4_${SAMPLE}.SbyNAme.bam
bamToFastq -i V4/V4_${SAMPLE}.SbyNAme.bam -fq V4/V4_${SAMPLE}_R1.fastq -fq2 V4/V4_${SAMPLE}_R2.fastq
gzip V4/V4_${SAMPLE}_R1.fastq
gzip V4/V4_${SAMPLE}_R2.fastq
rm V4/V4_${SAMPLE}.SbyNAme.bam
rm V4/V4_${SAMPLE}.sorted.bam

#############
echo "
merging V4 and unmapped datasets ..."
#############
mkdir V4_and_unmapped

cat V4/V4_${SAMPLE}_R1.fastq.gz unmapped/unmapped.${SAMPLE}_R1.fastq.gz > V4_and_unmapped/V4all_${SAMPLE}_R1.fastq.gz

cat V4/V4_${SAMPLE}_R2.fastq.gz unmapped/unmapped.${SAMPLE}_R2.fastq.gz > V4_and_unmapped/V4all_${SAMPLE}_R2.fastq.gz

#############
echo "
extracting reads that mapped PleD ..."
# The extracted reads are only those where both paired-end mapped
#############

samtools view -h -f 2 -F 4 ${SAMPLE}.bam Ga0074115_106:112275-112783 > PleD/PleD_${SAMPLE}.sorted.bam
samtools sort -n PleD/PleD_${SAMPLE}.sorted.bam -o PleD/PleD_${SAMPLE}.SbyNAme.bam
bamToFastq -i PleD/PleD_${SAMPLE}.SbyNAme.bam -fq PleD/PleD_${SAMPLE}_R1.fastq -fq2 PleD/PleD_${SAMPLE}_R2.fastq
gzip PleD/PleD_${SAMPLE}_R1.fastq
gzip PleD/PleD_${SAMPLE}_R2.fastq
rm PleD/PleD_${SAMPLE}.SbyNAme.bam
rm PleD/PleD_${SAMPLE}.sorted.bam

#############
echo "
extracting reads that mapped LpxA ..."
# The extracted reads are only those where both paired-end mapped
#############

samtools view -h -f 2 -F 4 ${SAMPLE}.bam Ga0074115_101:121173-121528 > LpxA/LpxA_${SAMPLE}.sorted.bam
samtools sort -n LpxA/LpxA_${SAMPLE}.sorted.bam -o LpxA/LpxA_${SAMPLE}.SbyNAme.bam
bamToFastq -i LpxA/LpxA_${SAMPLE}.SbyNAme.bam -fq LpxA/LpxA_${SAMPLE}_R1.fastq -fq2 LpxA/LpxA_${SAMPLE}_R2.fastq
gzip LpxA/LpxA_${SAMPLE}_R1.fastq
gzip LpxA/LpxA_${SAMPLE}_R2.fastq
rm LpxA/LpxA_${SAMPLE}.SbyNAme.bam
rm LpxA/LpxA_${SAMPLE}.sorted.bam

#############
echo "
extracting reads that mapped TufB ..."
# The extracted reads are only those where both paired-end mapped
#############

samtools view -h -f 2 -F 4 ${SAMPLE}.bam Ga0074115_122:16103-16402 > TufB/TufB_${SAMPLE}.sorted.bam
samtools sort -n TufB/TufB_${SAMPLE}.sorted.bam -o TufB/TufB_${SAMPLE}.SbyNAme.bam
bamToFastq -i TufB/TufB_${SAMPLE}.SbyNAme.bam -fq TufB/TufB_${SAMPLE}_R1.fastq -fq2 TufB/TufB_${SAMPLE}_R2.fastq
gzip TufB/TufB_${SAMPLE}_R1.fastq
gzip TufB/TufB_${SAMPLE}_R2.fastq
rm TufB/TufB_${SAMPLE}.SbyNAme.bam
rm TufB/TufB_${SAMPLE}.sorted.bam

done


##########
echo " 
Merging and aligning reads...
Here, bbtools is needed (https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/installation-guide/)"

##########
for file in `ls *.bam`
do
SAMPLE="${file/%.bam/}"


declare -a arr=("PleD" "LpxA" "TufB")
for locus in "${arr[@]}"
do

bash bbmerge.sh  qtrim2=r minq=28 in=${locus}/${locus}_${SAMPLE}_R1.fastq.gz in2=${locus}/${locus}_${SAMPLE}_R2.fastq.gz out=${locus}/${locus}_${SAMPLE}_merged.fastq.gz

##########
echo " 
Aligning reads with muscle...
For this part, the dependencies are: Seqtk (https://github.com/lh3/seqtk) and muscle (https://www.drive5.com/muscle/manual/install.html) "
##########

seqtk seq -a ${locus}/${locus}_${SAMPLE}_merged.fastq.gz > ${locus}/${locus}_${SAMPLE}_merged.fa

cat ../${locus}.fasta ${locus}/${locus}_${SAMPLE}_merged.fa > ${locus}/${locus}_${SAMPLE}_merged.fasta

rm ${locus}/${locus}_${SAMPLE}_merged.fa
muscle -in ${locus}/${locus}_${SAMPLE}_merged.fasta -out ${locus}/${locus}_${SAMPLE}_merged.aln

done



done
rm *.bai

echo "
Finished!"


cd ../