import CRISPR_functions
import pandas as pd
import numpy as np
import sys
import time
import Bio
import os
import argparse
import multiprocessing

if __name__=='__main__':
    
    parser = argparse.ArgumentParser(
        description = '** CRISPR_distance ** This script will take a list of CRISPR arrays and apply the probabilistic ordered independent loss model of Kupczok and Bollback (2013). The input file should have each array on a separate line. Each array should be represented by a comma separated list of spacer ids. \n\n\n \n DISCLAIMER: I am not a bioinformatician by any means so I may not be equipped to solve some of the issues that will arise but if you write to me, I will do by best. \n Send your comments and questions at https://github.com/maepz/CRISPR_distance or maeva.perez@umontreal.ca.',
        usage = 'Run_Parallel_CRISPR_Distance.py [options] <crispr arrays> <output directory>')
    parser.add_argument('inputfile', help='CRISPR arrays list input file', metavar='<inputfile>')
    parser.add_argument('outdir', help='name of the output directory', metavar='<outdir>')
    parser.add_argument("-n", help="number of processors to use. If not provided will use all available",dest='nproc')
    options = parser.parse_args()
    
    filename = options.inputfile
    dir = options.outdir
    if options.nproc:
        nproc = int(options.nproc)
    else:
        nproc = multiprocessing.cpu_count()
    print(f'number of processeors: {nproc}')
    
    if not os.path.isdir(dir):
        os.mkdir(dir)
    with open(filename,'r') as f:
        lines = f.read().splitlines()
    arrays = [line.split(',') for line in lines]
    start=time.time()
    Tree,dist,labels,LOG,final_dist = CRISPR_functions.main(arrays=arrays,nproc=nproc)
    labels_df = pd.DataFrame.from_dict(labels, orient="index")

    
    Bio.Phylo.write(Tree,dir+'/'+filename.split('/')[-1].replace('.txt','.nwk'),'newick')
    dist.to_csv(dir+'/'+filename.split('/')[-1].replace('.txt','_pairwiseDist.csv'),sep='\t')
    labels_df.to_csv(dir+'/'+filename.split('/')[-1].replace('.txt','_labels.csv'),sep='\t')
    LOG.to_csv(dir+'/'+filename.split('/')[-1].replace('.txt','_LogLikelihoods.csv'),sep='\t')
    
    # Export distances t1t2 to ancestor
    
    names=[i for i in range(len(arrays))]
    seqs=['-'.join(map(str,arr)) for arr in arrays]

    dic=dict(zip(names+seqs,seqs+names))
    dic
    df=pd.DataFrame(np.nan,columns=names,index=names)

    for k,v in final_dist.items():
        x=dic['-'.join(map(str,k.s1))] #row "distance from common ancestor to x..."
        y=dic['-'.join(map(str,k.s2))] # column "...when associated to y"
        d_xy=v[0]
        d_yx=v[1]

        df.loc[x,y]=d_xy
        df.loc[y,x]=d_yx
    df.to_csv(dir+'/'+filename.split('/')[-1].replace('.txt','_ancestorsDist.csv'),sep='\t')
    
    print('time to completion:',time.time()-start, 'seconds')
