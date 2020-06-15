import CRISPR_functions
import pandas as pd
import sys
import time

if __name__=='__main__':
    
    filename = sys.argv[1]

    with open(filename,'r') as f:
        lines = f.read().splitlines()
    arrays = [line.split(',') for line in lines]
    start=time.time()
    Tree,dist,labels = CRISPR_functions.main(arrays=arrays)
    dist.to_csv('Output/'+filename.split('/')[-1].replace('.txt','.csv'),sep='\t')
    labels_df = pd.DataFrame.from_dict(labels, orient="index")
    labels_df.to_csv('Output/'+'labels.csv')
    
    print('time to completion:',time.time()-start, 'seconds')