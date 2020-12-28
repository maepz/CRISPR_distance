#! /bin/bash

## Needs python 3+, and the following modules installed: argparse, re, os, sys, Biopython,  fuzzysearch, collections, pandas, numpy, gzip

echo"
## 1. download data and extract reads
"

python download_Amplicons_data.py -dir ./Amplicons_fastq

bash extract_merge_and_align_Amplicons_reads.sh ./Amplicons_fastq
echo"
## 2. detect haplotypes in TufB, PleD, and LpxA read alignments
# 2.1 Find SNPs. This step requires bowtie2 (http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml), samtools (http://www.htslib.org/doc/samtools.html), and VarScan (provided)
"

bash call_Amplicons_SNPs.sh ./Amplicons_fastq

echo"
# 2.2 Call haplotypes on detected SNP positions
"

python call_Amplicons_haplotypes.py -i Amplicons_fastq -o Amplicons_haplos
echo"
## 3. run dada2 script for V4 haplotype detection
Run the script dada2_V4.R in Rstudio for the haplotype detection based on the 16S V4 marker...
Press any key to continue"
while [ true ] ; do
read -t 3 -n 1
if [ $? = 0 ] ; then
exit ;
fi
done

