
import CRISPR_functions
from scipy.optimize import minimize,Bounds
from CRISPR_functions import neg_LL_floating_t
from CRISPR_functions import OPTIMIZE_pair_t1t2,OPTIMIZE_rho
from CRISPR_functions import is_overlapping,get_limits_ancestor_sizes,CRISPR_pair
from CRISPR_functions import get_distance_matrix_from_phylogeny
from CRISPR_functions import phylogeny_from_CRISPR
from Bio import Phylo
from Bio.Phylo import BaseTree,draw_ascii
import os,sys,time,itertools,multiprocessing
import numpy as np
import pandas as pd
from mpmath import mp,mpf,findroot

def get_nproc():
    # if on beluga, use environment variable cpu count
    if os.getenv('SLURM_CPUS_PER_TASK'):
        return(int(os.getenv('SLURM_CPUS_PER_TASK')))
    else:
        return(multiprocessing.cpu_count())
    
def get_ancestor_dict(args):
    '''Get dictionary of possible ancestors.'''
    pair,size_lims = args
    PAIR=CRISPR_pair(pair[0],pair[1])
    return(PAIR.get_combi(size_lims))

def main(arrays):

    # get number of cores
    nproc = get_nproc()
    #set precision
    mp.dps = 100; mp.pretty = True
    # initialize rho
    rho_init=mpf(sum([len(ar) for ar in arrays])/len([len(ar) for ar in arrays]))
    # overlapping pairs (list)
    overlapping_arrays=[pair for pair in list(itertools.combinations(arrays,2))
                        if is_overlapping(pair[0],pair[1])==1]
    # non-overlapping pairs (list)
    non_overlapping_arrays=[pair for pair in list(itertools.combinations(arrays,2))
                            if is_overlapping(pair[0],pair[1])==0]
    # pair list (list of PAIR objects)
    pairList = [CRISPR_pair(pair[0],pair[1]) for pair in overlapping_arrays]
    # size limits
    size_lims=get_limits_ancestor_sizes(arrays)
    # generate ancestor dictionary
    argList = [tuple([array,size_lims]) for array in overlapping_arrays]
    with multiprocessing.Pool(processes=nproc) as pool:
        ancestor_dicts = pool.map(get_ancestor_dict,argList)

    # first iteration: get initial t1s,t2s, and rho
    with multiprocessing.Pool(processes=nproc) as pool:
        # create argList for optimization of t1 and t2
        argList = [(pair,rho_init,size_lims,ancestor_dicts[index])
                   for index,pair in enumerate(pairList)]
        t1t2_list = pool.map(OPTIMIZE_pair_t1t2,argList)
    rho_update, neg_LL_update = OPTIMIZE_rho(t1t2_list,pairList,size_lims,non_overlapping_arrays,ancestor_dicts)

    Init_LL=[np.inf]
    i=0
    rho_list=[rho_update]
    convergence='False'
    
    # optimization of rho, and (t1, t2) of each pair
    while convergence=='False':
        i+=1
        print('iteration '+str(i)+'...')
        previous_LL=Init_LL[-1]
        # optimize t1 and t2
        with multiprocessing.Pool(processes=nproc) as pool:
            argList = [(pair,rho_update,size_lims,ancestor_dicts[index])
                       for index,pair in enumerate(pairList)]
            t1t2_list = pool.map(OPTIMIZE_pair_t1t2,argList)
        # optimize rho and return neg log likelihood
        rho_update, neg_LL_update = OPTIMIZE_rho(t1t2_list,pairList,size_lims,non_overlapping_arrays,ancestor_dicts)
        # check for convergence
        print(Init_LL, neg_LL_update, rho_update, rho_list)
        Init_LL+=neg_LL_update[0]
        rho_list+=rho_update
        if neg_LL_update[0]>previous_LL:
            convergence='True'
    final_dist=dict(zip(pairList,t1t2_list))
    
    # Get phylogeny and distance matrix
    Tree,labels=phylogeny_from_CRISPR(arrays,final_dist)
    dist=get_distance_matrix_from_phylogeny(Tree)
    print(f'Final rho: {rho_list[-1]}')
    return(Tree,dist,labels)
    
    
if __name__=='__main__':
    
    filename = sys.argv[1]
    
    with open(filename,'r') as f:
        lines = f.read().splitlines()
    arrays = [line.split(',') for line in lines]

    Tree,dist,labels = main(arrays=arrays)
