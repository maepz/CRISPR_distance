#!/usr/bin/env python
# coding: utf-8

# In[2]:


import argparse
import re
from Bio import Entrez
import os
import sys


# In[ ]:


if __name__=='__main__':
    
    parser = argparse.ArgumentParser(
        description = '** Download CRISPR data ** This script will download the CRISPR PacBio reads from the Bioproject PRJNA641184 as fasta.gz files into a chosen directory ',
        usage = 'download_CRISPR_data.py -dir <directory>')
    parser.add_argument('-dir', help='name of the directory in which the reads will be downloaded; if it doesnt exist, it will be created', dest='dirname')
    options = parser.parse_args()
    
    dirname = options.dirname
    if not os.path.isdir(dirname):
        os.mkdir(dirname)

## Fetch the SRAs from the Bioproject
Entrez.email = "maeva.perez@umontreal.ca"  # Always tell NCBI who you are
handle = Entrez.esearch(db="sra", term="PRJNA641184[BioProject] AND PacBio_SMRT[PLATFORM]", retmax=100)
search_results = Entrez.read(handle)
print('There are %s CRISPR SRAs' % len(search_results['IdList']))
print('SRA ids: ',search_results['IdList'])

## Download the PacBio raw reads for each SRA; necessitates SRA tools installed

i=0
starturl='url='
endurl='size'
startlib='<LIBRARY_NAME>'
endlib='</LIBRARY_NAME>'
srrs=[]
for result in search_results["IdList"]:
    handle = Entrez.efetch(db="sra", id=result, retmode='xml')
    srr = handle.read()
    url = re.search('%s(.*)%s' % (starturl, endurl), srr).group(1)[1:-2]
    lib = re.search('%s(.*)%s' % (startlib, endlib), srr).group(1)
    srr_id = url.split('/')[-1]
    srrs+=[srr_id]
    print( "Downloading "+ srr_id + "; sample " + lib + "...")
    cmd = "prefetch "+ srr_id
    os.system(cmd)
    print( "\nDownload Complete "+ srr_id)
    
    
    print( "\nConverting...")   
    cmd = 'fasterq-dump '+srr_id+' ; gzip '+srr_id+'.fastq ; mv '+srr_id+'.fastq.gz '+dirname+'/'+lib+'.fastq.gz'
    os.system(cmd)
    
    cmd = 'rm -r '+srr_id
    os.system(cmd)
    
    i+=1
#    if i> 1:
#        break
