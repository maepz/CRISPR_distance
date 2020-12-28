#!/usr/bin/env python
# coding: utf-8

## Load necessary packages
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import Gapped,IUPAC
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import generic_dna
from os import listdir
from os.path import isfile, join
import pandas as pd
from collections import Counter
import argparse
import os

if __name__=='__main__':

    parser = argparse.ArgumentParser(
        description = '** call_Amplicons_haplotypes ** \nThis script will call Amplicons haplotypes based on genotype at predefined variable positions. Output files are written to a chosen direcory.',
        usage = 'call_Amplicons_haplotypes.py -i <input directory> -o <output directory>')
    parser.add_argument('-i', help='input directory containing the subdirectories for each amplicons', dest='inputdir')
    parser.add_argument('-o', help='name of the output directory for the output files; if it doesnt exist, it will be created', dest='outdir')

    options = parser.parse_args()

    inputdir = options.inputdir+'/'
    outdir = options.outdir+'/'

    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    
    
########## FIND HAPLOTYPES IN LpxA ##############
print('Calling haplotypes for LpxA...')

LpxA_df=pd.DataFrame(columns=['Sample'])
seqlist=[]

mypath=inputdir+'LpxA/'
onlyfiles = [mypath+f for f in listdir(mypath) if f[-10:]=='merged.aln']

for file in onlyfiles:

    sample=file.split('/')[-1][5:-11]
    align = AlignIO.read(file, "fasta",alphabet=Gapped(IUPAC.ambiguous_dna))
    aln = MultipleSeqAlignment([f for f in sorted(align, key=lambda x : x.id)])
    

    #### remove ref gaps in alignment ####
    gap=aln[0].seq.find('-')
    while gap>-1:
        aln=aln[:, :gap] + aln[:, gap+1:]
        gap=aln[0].seq.find('-')
    snp1=aln[0].seq.find('CTTGTTGCTGACGG')-1
    snp2=aln[0].seq.find('TCTCGATCACGCAG')-1
    dic=Counter([''.join(list(i)) for i in zip(aln[:, snp1],aln[:, snp2])])
    dic['Sample']=sample
    LpxA_df=LpxA_df.append(dic, ignore_index=True)

# remove artefactual haplotypes
for col in list(LpxA_df)[1:]:
    if LpxA_df[col].sum()<len(LpxA_df['Sample']):
        LpxA_df=LpxA_df.drop(col,axis=1)

LpxA_df=LpxA_df.add_prefix('LpxA_')
LpxA_df.rename(columns={'LpxA_Sample':'Sample'},inplace=True)

LpxA_df.to_csv(outdir+'LpxA_haplotypes.txt',header=True,sep='\t',index=False)



########## FIND HAPLOTYPES IN PleD ##############
print('Calling haplotypes for PleD..')

PleD_df=pd.DataFrame(columns=['Sample'])
seqlist=[]


mypath=inputdir+'PleD/'
onlyfiles = [mypath+f for f in listdir(mypath) if f[-10:]=='merged.aln']
for file in onlyfiles:

    sample=file.split('/')[-1][5:-11]
    align = AlignIO.read(file, "fasta",alphabet=Gapped(IUPAC.ambiguous_dna))
    aln = MultipleSeqAlignment([f for f in sorted(align, key=lambda x : x.id)])
    

    #### remove ref gaps in alignment ####
    gap=aln[0].seq.find('-')
    while gap>-1:
        aln=aln[:, :gap] + aln[:, gap+1:]
        gap=aln[0].seq.find('-')
    snp1=aln[0].seq.find('GTATTCTGCTCACC')-1
    snp2=aln[0].seq.find('GGCCGGCAGATGCA')-1
    dic=Counter([''.join(list(i)) for i in zip(aln[:, snp1],aln[:, snp2])])
    dic['Sample']=sample
    PleD_df=PleD_df.append(dic, ignore_index=True)

for col in list(PleD_df)[1:]:
    if PleD_df[col].sum()<len(PleD_df['Sample']):
        PleD_df=PleD_df.drop(col,axis=1)

        
PleD_df = PleD_df.add_prefix('PleD_')
PleD_df.rename(columns={'PleD_Sample':'Sample'},inplace=True)

PleD_df.to_csv(outdir+'PleD_haplotypes.txt',header=True,sep='\t',index=False)


########## FIND HAPLOTYPES IN TufB ##############
print('Calling haplotypes for TufB..')

TufB_df=pd.DataFrame(columns=['Sample'])
seqlist=[]


mypath=inputdir+'TufB/'
onlyfiles = [mypath+f for f in listdir(mypath) if f[-10:]=='merged.aln']
for file in onlyfiles:

    sample=file.split('/')[-1][5:-11]
    align = AlignIO.read(file, "fasta",alphabet=Gapped(IUPAC.ambiguous_dna))
    aln = MultipleSeqAlignment([f for f in sorted(align, key=lambda x : x.id)])
    

    #### remove ref gaps in alignment ####
    gap=aln[0].seq.find('-')
    while gap>-1:
        aln=aln[:, :gap] + aln[:, gap+1:]
        gap=aln[0].seq.find('-')
    snp1=aln[0].seq.find('ATGATCACAGGTGT')-1
    snp2=aln[0].seq.find('GGCAGAGACCACCA')-1
    dic=Counter([''.join(list(i)) for i in zip(aln[:, snp1],aln[:, snp2])])
    dic['Sample']=sample
    TufB_df=TufB_df.append(dic, ignore_index=True)

for col in list(TufB_df)[1:]:
    if TufB_df[col].sum()<len(TufB_df['Sample']):
        TufB_df=TufB_df.drop(col,axis=1)
        
PleD_df = PleD_df.add_prefix('TufB_')
PleD_df.rename(columns={'TufB_Sample':'Sample'},inplace=True)

TufB_df.to_csv(outdir+'TufB_haplotypes.txt',header=True,sep='\t',index=False)


### All data
df=pd.merge(LpxA_df, PleD_df, on='Sample')
df=pd.merge(df, TufB_df, on='Sample')
df.to_csv(outdir+'haplotypes.txt',header=True,sep='\t',index=False)

