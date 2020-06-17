import CRISPR_functions
import pandas as pd
import numpy as np
import sys
import time
import Bio

if __name__=='__main__':
    
    filename = sys.argv[1]

    with open(filename,'r') as f:
        lines = f.read().splitlines()
    arrays = [line.split(',') for line in lines]
    start=time.time()
    Tree,dist,labels,LOG,final_dist = CRISPR_functions.main(arrays=arrays)
    labels_df = pd.DataFrame.from_dict(labels, orient="index")

    
    Bio.Phylo.write(Tree,'Output2/'+filename.split('/')[-1].replace('.txt','.nwk'),'newick')
    dist.to_csv('Output2/'+filename.split('/')[-1].replace('.txt','_pairwiseDist.csv'),sep='\t')
    labels_df.to_csv('Output2/'+filename.split('/')[-1].replace('.txt','_labels.csv'),sep='\t')
    LOG.to_csv('Output2/'+filename.split('/')[-1].replace('.txt','_LogLikelihoods.csv'),sep='\t')
    
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
    df.to_csv('Output2/'+filename.split('/')[-1].replace('.txt','_ancestorsDist.csv'),sep='\t')
    
    print('time to completion:',time.time()-start, 'seconds')
