#!/usr/bin/env python
# coding: utf-8

# In[10]:

from fuzzysearch import find_near_matches

import os
from os import listdir
from os.path import isfile, join

from Bio import SeqIO
from Bio import AlignIO
from Bio import SeqRecord
from Bio import Seq
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from collections import Counter

import pandas as pd
import numpy as np
import sys
import gzip
import argparse



# In[ ]:


if __name__=='__main__':
    
    parser = argparse.ArgumentParser(
        description = '** call_CRISPR_haplotypes ** \nThis script will call CRISPR haplotypes based on near-matching the known spacers. Version 6. Both primers have to be present in the reads for haplotype calling. Output files are written to a chosen direcory. The outputs are: \nALLCRISPRs_haplo_matches_at_5_dist_V6_sorted.txt = A csv file with the haplotype counts in all samples;\n SAMPLE.spacers_matches_v6.txt = A table with the matching pattern for each read in SAMPLE;\n SAMPLE.Good_records.fasta = A fasta file containing all reads assigned to a CRISPR haplotype;\n SAMPLE.Missing_spacers_records.fasta = A fasta file containing all other reads',
        usage = 'call_CRISPR_haplotypes.py -i <input directory> -o <output directory>')
    parser.add_argument('-i', help='input directory containing the fastq.gz files', dest='inputdir')
    parser.add_argument('-o', help='name of the output directory for the output files; if it doesnt exist, it will be created', dest='outdir')
    
    options = parser.parse_args()
    
    mypath = options.inputdir+'/'    
    outdir = options.outdir+'/'
    
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

# Define spacers to find 

flanker=['flanker','TCCGATCGGCCAAAAAGATCGGTAGAACTCAACACGGTGATTTTTCTTATCTCTATCAAATAGATAGAATAAGA']
flk=['flk',flanker[1][-32:]]
sp5=['sp5','GGCTGGTGGTGCTGACGATGGCGATCATCCAC']
sp170=['sp170','AACATTAACTGCTTGCCTTTCGGCTTTTTTGT']
sp51=['sp51','TGTTTGGCATGGACTGGATCGCGAGCGGTGG']
sp22=['sp22','GGACCGTGATCCTGACCGGCGCTGCGCTTGA']
sp79=['sp79','CTCTACAACACGCCACTTCTTTGTCGTTTGA']
sp161=['sp161','CCACCCCGGAGATAGTTACCGATTCCTGGCGG']
sp174=['sp174','GCAATTGCCAAGATTGGCTTTGCCCTGGTACT']
sp179=['sp179','TGGTCGTTTCGTGCTGGCTACGGGGTTGGGTC']
sp37=['sp37','ATTCCGGGCCGCAAATATCCGACGGCTTAGAA']
sp47=['sp47','CATTCTTTCCGTTCTGGCGGATGATCCATGTG']
sp440=['sp440','GTCAATATTCCGCGCCCACAAAGGCGTTATTG']
sp184=['sp184','GCACTCGGCCAGTCGACGGACATAGGTGCGAT']
sp55=['sp55','CATCAATGAGCGCCAGTGTTCTGACGCCGTAA']
sp58=['sp58','TTTGGACTGCCTCGTATCCCTTGATAATCGCC']
sp295=['sp295','ACTTACGATCTGATTGCCTCGCGTAATAACT']
tailprimer=['tailprimer','TCTTGAAGCGAGGTCACTCTCCGGGTACAGGC']

spacers=[flk,sp5,sp170,sp51,sp22,sp79,sp161,sp174,sp179,sp37,sp47,sp440,sp184,sp55,sp58,sp295,tailprimer]
        

# init haplo_table

Haplo_matches=pd.DataFrame(0, index=[],columns=[])
output_table ='ALLCRISPRs_haplo_matches_at_5_dist_V6_sorted.txt'

# Spacer match algorithm V6

for file in listdir(mypath):
    with gzip.open(mypath+file, "rt") as handle:
        sequences=[seq for seq in SeqIO.parse(handle, "fastq")]
        
    sample_id=file[4:-9]
    
    print( 'Working on '+sample_id+'...')
    
    #### init table matches and other params ####
    idx=[seq.id for seq in sequences]
    cols=[sp[0] for sp in spacers]+['seqlen']+['haplo']+['warnings']+['sep_lengths']
    Seq_matches=pd.DataFrame(str(), index=idx, columns=cols)
    
    dic={}
    maxdist=5 #maximum (hamming) distance between spacer and read
    totrec=0
    rec=0
    Matched_records=[]
    Missing_spacers_records=[]
    

    #### read fastq and fill tables ####
#     sequences=(seq for seq in SeqIO.parse(mypath+file+'.fastq', "fastq"))
    for seq in sequences:
        totrec+=1
        dic[seq.id]=[]
        Seq_matches.loc[seq.id]['seqlen']=len(seq)
        
        
        ### ignore sequences that are too small
        if len(seq.seq)<33:
            Seq_matches.loc[seq.id]['warnings']=['blast_not_run']
            Seq_matches.loc[seq.id]['haplo']='No_matches'
            continue # go to next read

        
        ### detect spacers ###
        for n in range(len(spacers)): # match each spacer sequence to read
            rc=''
            match=find_near_matches(spacers[n][1], str(seq.seq), max_l_dist=maxdist) #match spacer to forward sequence
            if len(match)>0: # if match, record the match
                dic[seq.id]+=[spacers[n][0],rc]
                Seq_matches.loc[seq.id][spacers[n][0]]=match
            if len(match)==0: #if no match found, try matching spacer to reverse-complement sequence
                match=find_near_matches(spacers[n][1], str(seq.seq.reverse_complement()), max_l_dist=maxdist)
                if len(match)>0:
                    rc='rc'
                    dic[seq.id]+=[spacers[n][0],rc]
                    Seq_matches.loc[seq.id][spacers[n][0]]=match                        


        ## check if all matched spacers are in the same orientation ##
        if 'rc' in dic[seq.id] and dic[seq.id][1::2] != len(dic[seq.id][::2])*['rc']:
            Seq_matches.loc[seq.id]['warnings']='not_all_spacers_are_in_the_same_orientation'
            pass
        elif 'rc' in dic[seq.id] and dic[seq.id][1::2] == len(dic[seq.id][::2])*['rc']:
            rcseq=seq.reverse_complement()
            rcseq.id=seq.id
            rcseq.description='rc'
            seq=rcseq
        
        
        ## add haplotype info to reads ##
        if len([s for s in Seq_matches.loc[seq.id].values if s!='']) == 0:
            Seq_matches.loc[seq.id]['haplo']='No_matches'
            continue
            
        else:
#             print([Seq_matches.loc[seq.id].index[i] for i in range(len(Seq_matches.loc[seq.id].index)-4)])
            Seq_matches.loc[seq.id]['haplo']='_'.join([Seq_matches.loc[seq.id].index[i] if len(Seq_matches.loc[seq.id].values[i])>0 else 'NA' for i in range(len(Seq_matches.loc[seq.id].index)-4)])

            if not all(spacers in dic[seq.id]  for spacers in ['flk','tailprimer']): #read missing flanker of end
                Seq_matches.loc[seq.id]['warnings']=list(Seq_matches.loc[seq.id]['warnings'])+['missing_anchors']
        
        ## check if matched spacers are separated by a single repeat ##
        row=[item[0] for item in Seq_matches.loc[seq.id,[u'flk', u'sp5', u'sp170', u'sp51', u'sp22', u'sp79',
       u'sp161', u'sp174', u'sp179', u'sp37', u'sp47', u'sp440', u'sp184',
       u'sp55', u'sp58', u'sp295']].dropna().values if item != '']
        try:
            sep_lengths = [np.abs(row[i+1].start-row[i].end) for i in range(len(row)-1)]
            maxx=max(sep_lengths)
            minn=min(sep_lengths)

        except ValueError:
            sep_lengths = [np.abs(row[i+1].start-row[i].end) for i in range(len(row)-1)]
            maxx=1
            minn=1
        maxrelmin=maxx/minn
        if maxrelmin>=2:
            Seq_matches.loc[seq.id]['warnings']=list(Seq_matches.loc[seq.id]['warnings'])+['some spacers missing']
#             Seq_matches.loc[seq.id]['haplo']='Imperfect_match'
            Seq_matches.loc[seq.id]['sep_lengths']=sep_lengths

            
        ## if read/record has no warnings, save it for outputting in fasta ##
        if len(Seq_matches.loc[seq.id]['warnings'])==0:
            rec+=1
            Matched_records+=[seq]
        if 'missing_anchors' not in Seq_matches.loc[seq.id,'warnings'] and 'some spacers missing' in Seq_matches.loc[seq.id,'warnings']:
            Missing_spacers_records+=[seq]

    ### write fasta, and fastq for matched records ###
    SeqIO.write(Matched_records, outdir+sample_id+'.Good_records.fasta', "fasta")
    SeqIO.write(Matched_records, outdir+sample_id+'.Good_records.fastq', "fastq")
    SeqIO.write(Missing_spacers_records, outdir+sample_id+'.Missing_spacers_records.fasta', "fasta")

    ### write summary table for sample ###
    Seq_matches.to_csv(path_or_buf=outdir+sample_id+'.spacers_matches_v6.txt', sep='\t', na_rep='', header=True, index=True, index_label='seq.id', mode='w', encoding=None, compression=None, quoting=None)


    ### fill haplo table ###
    ## parse Seq_matches to keep only those with no warnings (those that have both anchors and all spacers with the same orientation and spacing) ##
    Seq_matches['warnings'].astype('str')
    mask= Seq_matches['warnings'].str.len() == 0 
    parsed_Seq_matches=Seq_matches[mask]
    counts=Counter(parsed_Seq_matches['haplo'])
    df = pd.DataFrame.from_dict(counts, orient='index',columns=[sample_id])
    Haplo_matches = pd.concat([Haplo_matches, df], axis=1, sort=True)

    ### print verbose ###
    print( 'finished parsing '+ sample_id+'\n'+    'total number of reads='+ str(totrec)+'\n'+    'number of reads with haplotypes='+str(rec)+' ('+str(100*rec/totrec)+'%)')
    
    

# Sort df to place main haplotypes on the top and write output

main_haplos = {0:'flk_sp5_sp170_sp51_sp22_sp79_sp161_sp174_sp179_sp37_sp47_sp440_sp184_sp55_sp58_sp295_tailprimer',
1:'flk_sp5_sp170_sp51_sp22_sp79_sp161_sp174_sp179_NA_NA_sp440_sp184_sp55_sp58_sp295_tailprimer',
7:'flk_sp5_sp170_NA_NA_NA_sp161_sp174_sp179_sp37_sp47_sp440_sp184_sp55_NA_sp295_tailprimer',
3:'flk_sp5_sp170_NA_NA_NA_sp161_sp174_sp179_sp37_NA_sp440_sp184_sp55_sp58_sp295_tailprimer',
4:'flk_sp5_sp170_NA_NA_NA_sp161_sp174_sp179_sp37_sp47_NA_sp184_sp55_sp58_sp295_tailprimer',
5:'flk_sp5_sp170_NA_NA_NA_sp161_sp174_sp179_sp37_sp47_sp440_sp184_NA_sp58_sp295_tailprimer',
6:'flk_sp5_sp170_NA_NA_NA_NA_sp174_sp179_sp37_sp47_sp440_sp184_sp55_NA_sp295_tailprimer',
2:'flk_sp5_sp170_NA_NA_NA_sp161_sp174_sp179_sp37_sp47_sp440_sp184_sp55_sp58_sp295_tailprimer',
8:'flk_sp5_sp170_NA_NA_NA_sp161_sp174_sp179_NA_sp47_sp440_sp184_sp55_NA_sp295_tailprimer',
9:'flk_sp5_sp170_NA_NA_NA_sp161_NA_NA_NA_NA_NA_sp184_sp55_NA_sp295_tailprimer'}

main_haplos = dict([[v,k] for k,v in main_haplos.items()])
main_haplos

Haplo_matches['new'] = Haplo_matches.index.map(main_haplos)
Haplo_matches = Haplo_matches.sort_values("new")
Haplo_matches = Haplo_matches.reset_index()
Haplo_matches = Haplo_matches.drop('new',axis=1)

Haplo_matches.to_csv(outdir+output_table)
print( 'All done!')
