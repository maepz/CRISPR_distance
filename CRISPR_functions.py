#!/usr/bin/env python
# coding: utf-8

# dependencies
import pandas as pd
from Bio.Phylo import BaseTree,draw_ascii
from Bio.Phylo.BaseTree import Clade,Tree
from scipy.optimize import minimize,Bounds
from scipy.stats import poisson
from mpmath import mpf,mp,findroot
import itertools,os,multiprocessing
import math
import numpy as np

# Here are the functions needed for CRISPR_distance

# ## Get the ansestors
def get_limits_ancestor_sizes(arrays,p_low=0.005,p_high=0.995):
    '''Given a set of arrays and rho, returns l1 and l2 the min and max ancestor lengths'''
    rho=sum([len(ar) for ar in arrays])/len([len(ar) for ar in arrays])
    return(tuple(poisson.ppf([p_low,p_high],rho).astype(int)))

def nCr(n,r):
        f = math.factorial
        return f(n) / f(r) / f(n-r)

def categorize_spacers_for_ordered_model(s1,s2):
    ''' Compare two CRISPR arrays and categorize spacers into four categories: 
    c1_spacers: spacers in s1 present in the ancestor, 
    c2_spacers: spacers in s2 present in the ancestor,
    d1_spacers: spacers unique to s1
    d2_spacers: spacers unique to s2'''
    c_spacers=[]
    c1_spacers=0
    c2_spacers=0
    d1_spacers=0
    d2_spacers=0
    for i in range(len(s1)):
        sp=s1[i]
        if sp in s2:
            j=s2.index(sp)
            c1_spacers=s1[i:]
            c2_spacers=s2[j:]
            d1_spacers=s1[:i]
            d2_spacers=s2[:j]
            break
    return(c1_spacers,c2_spacers,d1_spacers,d2_spacers)

def is_overlapping(s1, s2):
    ''' Are two CRISPR arrays overlapping? '''
    for i in s1:
        if i in s2:
            answer=1
            break
        else:
            answer=0
            continue
    return(answer)


class combi:
    
    """The combi object is a list of integers [c,i,j,u] that represent all putative ancestors given the following parameters:
    c = number of c_spacers; spacers necessarily present in the ancestor of two CRISPR arrays,
    i = number of putative ancestral spacers amongst the spacers only present in array1,
    j = number of putative ancestral spacers amongt these only present in array 2,
    u = number of putative spacers present in ancestor but lost in both arrays
    
    the lengths of the putative ancestral arrays n = sum([c,i,j,u])"""
    
    def __init__(self,liste):
        self.list = liste
        self.c = liste[0]
        self.i = liste[1]
        self.j = liste[2]
        self.u = liste[3]
        self.n = sum(self.list)
        possible_positions_of_u = nCr(self.n,self.u)
        possible_combinations_of_d1_and_d2=nCr(self.i+self.j,self.i)
        
        self.array_counts = possible_positions_of_u*possible_combinations_of_d1_and_d2 # to get the number of putative ancestors
        self.str = '-'.join(map(str,self.list))
        
    def compare_ancestor_to_arrays(self,PAIR):
        '''combi_s0= combi class object, PAIR= CRISPR_pair class object, returns for each array
        m number of spacers maintained between ancestor and array
        d number of spacers lost in array
        j number of spacers gained in array'''

        c=self.c
        i=self.i
        j=self.j
        u=self.u
        c1=len(PAIR.c1_spacers)
        c2=len(PAIR.c2_spacers)
        dd1=len(PAIR.d1_spacers)
        dd2=len(PAIR.d2_spacers)
        
        m1=c1+i
        d1=j+u+c-c1
        j1=dd1-i
        m2=c2+j
        d2=i+u+c-c2
        j2=dd2-j

        return(m1,d1,j1,m2,d2,j2)
        
    def get_arrays(self,c_spacers,d1_spacers,d2_spacers): # to get the list of all possible arrays not usefull       
        array_list2=[]
        def merge_lists(lst1,lst2):
            
            array_list=[]
            for locations in itertools.combinations(range(len(lst1) + len(lst2)), len(lst2)):
                result = lst1[:]
                for location, element in zip(locations, lst2):
                    result.insert(location, element)
                new_list=result
                array_list+=[new_list]
            return(array_list)
        
        d1_good=d1_spacers[len(d1_spacers)-self.i:]
        d2_good=d2_spacers[len(d2_spacers)-self.j:]
        u_good=self.u*['u']
        c_spacers=c_spacers
        if self.u>0:
            for array in merge_lists(d1_good,d2_good):
                array_list2+=merge_lists(u_good,array+c_spacers)
        else:
            array_list2=[array+c_spacers for array in merge_lists(d1_good,d2_good)]
        return(array_list2)
    

class CRISPR_pair:
    
    def __init__(self, s1,s2):
        self.s1 = s1
        self.s2 = s2
        self.c1_spacers, self.c2_spacers, self.d1_spacers, self.d2_spacers = categorize_spacers_for_ordered_model(self.s1,self.s2)
        self.c_spacers = list(set(self.c1_spacers+self.c2_spacers))
        
        
    def get_combi(self,size_lims):

        ''' The function get_combi outputs a dictionary of all the possible combinations of spacers categories and their corresponding adjusted (per ancestor length) weights. 
        {c-i-j-u:ws} with:
        c number of spacers in common (spacers necessarily present in ancerstor), 
        i number of ancestral spacers amongst the spacers only present in array1, 
        j number of ancestral spacers amongt these only present in array2,
        u number of ancestral spacers lost in both array1 and array2,
        ws weight of each putative ancestral array from this combi.
        l1 (min ancestor length) and l2 (max ancestor length) have to be provided
        n length of ancestral array = sum(c,i,j,u) '''

        mp.dps = 100; mp.pretty = True
       
        l1=size_lims[0]
        l2=size_lims[1]
        min_ansestor_len=min(len(self.c_spacers),l1)
        max_ansestor_len=max(len(self.c_spacers+self.d1_spacers+self.d2_spacers),l2)
        c=len(self.c_spacers)
        spacers_combi={}
        combi_output={}
        for n in range(min_ansestor_len,max_ansestor_len+1):
            for i in range(len(self.d1_spacers)+1):
                for j in range(len(self.d2_spacers)+1):
                    if c+i+j<n+1:
                        u=n-(c+i+j)
                        if n in spacers_combi.keys():
                            spacers_combi[n]+=[combi([c,i,j,u])]
                        else:
                            spacers_combi[n]=[combi([c,i,j,u])]
        
        for n,v in spacers_combi.items():
            
            combi_strings=[combi.str for combi in v]
            ancestor_counts = sum([combi.array_counts for combi in v])
            ancestor_weights = [mpf(combi.array_counts)/ancestor_counts for combi in v]
            
            combi_output.update(dict(zip(combi_strings,ancestor_weights)))

        return(combi_output)


## Math equations

### Length model 

def rho(lbda,mu):
    return(lbda/mu)

def prob_n_given_ro(n,rho):
    return((np.e**-rho)*(rho**n/np.math.factorial(n))) 

### Transitions prob for the ordered model

def M(m,t,mu):
    '''probability of preserving m spacers in time t : M(m|t,mu)'''
    mp.dps = 100; mp.pretty = True
    exp_array = np.frompyfunc(mp.exp, 1, 1)
    M=exp_array(-m*t*mu)
    return(M)

def D(d,t,mu):
    '''probability of loosing d spacers in time t : D(d|t,mu)'''
    mp.dps = 100; mp.pretty = True
    exp_array = np.frompyfunc(mp.exp, 1, 1) # to make mp.exp work with arrays
    power_array = np.frompyfunc(mp.power, 2, 1) # to make mp.power work with arrays
    D=power_array(1 - exp_array(-t*mu),d)
    return(D)

def I(j,t,lbda,mu):
    '''probability of inserting d spacers in time t : D(d|t,mu)'''
    mp.dps = 100; mp.pretty = True
    rho=lbda/mu
    exp_array = np.frompyfunc(mp.exp, 1, 1)
    power_array = np.frompyfunc(mp.power, 2, 1)
    a=rho*(1 - exp_array(-t*mu))
    I=(exp_array(-a)*power_array(a,j))/np.math.factorial(j)
    return(I)


def L(rho,t,PAIR,size_lims,ancestor_dict=None):
    '''The likelihood of a pair of spacer arrays (s1, s2) with
times (t1, t2) anf rho'''
    t1=t[0]
    t2=t[1]
    l1=size_lims[0]
    l2=size_lims[1]
    lbda=rho/(rho+1)
    mu=1/(rho+1)

    Likelihood=0
    
    for s0,ws in ancestor_dict.items():
        combi_s0 = combi(list(map(int,s0.split('-'))))
        n=combi_s0.n

        m1,d1,j1,m2,d2,j2 = combi_s0.compare_ancestor_to_arrays(PAIR)

        Qa = prob_n_given_ro(n,rho)*ws
        T1 = M(m1,t1,mu)*D(d1,t1,mu)*I(j1,t1,lbda,mu)
        T2 = M(m2,t2,mu)*D(d2,t2,mu)*I(j2,t2,lbda,mu)
        Likelihood += Qa*T1*T2

        if Qa*T1*T2<0:
            print('rho',rho)

    return(Likelihood)

def neg_LL_floating_t(x,rho,PAIR,size_lims,ancestor_dict):
    mp.dps = 100; mp.pretty = True

    log_array = np.frompyfunc(mp.log, 1, 1) # to make mp.log work with arrays
    return(np.float(-log_array(L(rho=rho,t=x,PAIR=PAIR,size_lims=size_lims,ancestor_dict=ancestor_dict))))
    
def neg_LL_floating_rho(x,t1t2_list,pair_list,size_lims,non_overlapping_arrays,ancestor_dicts):
    mp.dps = 100; mp.pretty = True
    
    log_array = np.frompyfunc(mp.log, 1, 1)

    neg_LL_overlapping=sum([-log_array(L(rho=x,t=t1t2_list[i],PAIR=pair_list[i],size_lims=size_lims,ancestor_dict=ancestor_dicts[i])) for i in range(len(pair_list))])
    neg_LL_non_overlapping=sum([-log_array(prob_n_given_ro(len(pair[0]),x)*prob_n_given_ro(len(pair[1]),x)) for pair in non_overlapping_arrays])
    
    return(neg_LL_overlapping+neg_LL_non_overlapping)


## Optimizers


def OPTIMIZE_rho(t1t2_list,pair_list,size_lims,non_overlapping_arrays,ancestor_dicts):
    '''Provided a set of divergence times t1,t2 from ancestor to arrays for all array pairs, OPTIMIZE_rho finds the best rho'''
    x0=[2]
    bounds=Bounds(lb=0,ub=np.inf)
#     optimize=minimize(neg_LL_floating_rho,x0,bounds=bounds,method='Powell',args=(t1t2_list,pair_list,size_lims,non_overlapping_arrays)) # Powell is overshooting at some point and try a rh0<0; it seems bounds have not been implemented for this method
    optimize=minimize(neg_LL_floating_rho,x0,bounds=bounds,args=(t1t2_list,pair_list,size_lims,non_overlapping_arrays,ancestor_dicts)) # method=L-BFGS-B
    return(optimize.x,optimize.fun)


def OPTIMIZE_t1t2(overlapping_arrays, rho, size_lims):
    '''!!! DEPRECADTED: USE INSTEAD THE PARALLEL VERSION !!!Provided a set of overlapping arrays and rho, OPTIMIZE_t1t2 finds their best respective divergence times t1,t2 from ancestor to arrays. The output is an array of [t1,t2] of length len(overlapping_arrays)'''
    t1t2_list=[]

    for pair in overlapping_arrays:
        PAIR=CRISPR_pair(pair[0],pair[1])
        ancestor_dict = PAIR.get_combi(size_lims)
        x0=[1,1]
        bnds = ((0, None), (0, None))
#         optimize=minimize(neg_LL_floating_t,x0,bounds=bounds,method='Powell',args=(rho,PAIR,size_lims))
        optimize=minimize(neg_LL_floating_t,x0,bounds=bnds,args=(rho,PAIR,size_lims,ancestor_dict)) # method=L-BFGS-B
        t1t2_list+=[tuple(optimize.x)]
    
    return(t1t2_list)

def OPTIMIZE_pair_t1t2(args):
    '''Provided an individual PAIR object, rho, size limits, and the pair ancestor dictionary, OPTIMIZE_t1t2 finds their best respective divergence times t1,t2 from ancestor to arrays. The output is (t1,t2) for the pair.'''
    
    PAIR, rho, size_lims, ancestor_dict = args

    x0=[1,1]
    bnds = ((0, None), (0, None))
#     optimize=minimize(neg_LL_floating_t,x0,bounds=bounds,method='Powell',args=(rho,PAIR,size_lims))
    optimize=minimize(neg_LL_floating_t,x0,bounds=bnds,args=(rho,PAIR,size_lims,ancestor_dict)) # method=L-BFGS-B
    return(tuple(optimize.x))

def find_optimum_rho_and_distances_ordered_model(arrays):
    '''
    !!! DEPRECATED: USE INSTREAD THE PARALLEL VERSION !!!
    Given a list of arrays as input, optimize rho and t for the ordered model
    
    Output:dictionnary of pairs and their respecive distances, rho
    '''
    overlapping_arrays=[pair for pair in list(itertools.combinations(arrays,2)) if is_overlapping(pair[0],pair[1])==1]

    ## first step
    rho_init=mpf(sum([len(ar) for ar in arrays])/len([len(ar) for ar in arrays]))
    size_lims=get_limits_ancestor_sizes(arrays)
    t1t2_list=OPTIMIZE_t1t2(overlapping_arrays, rho_init, size_lims)
    pair_list=[CRISPR_pair(pair[0],pair[1]) for pair in overlapping_arrays]
    ancestor_dicts = [PAIR.get_combi(size_lims) for PAIR in pair_list]
    non_overlapping_arrays=[pair for pair in list(itertools.combinations(arrays,2)) if is_overlapping(pair[0],pair[1])==0]
    (rho_update, neg_LL_update) = OPTIMIZE_rho(t1t2_list,pair_list,size_lims,non_overlapping_arrays,ancestor_dicts)
    ## iterate
    Init_LL=[np.inf]
    i=0
    step_list=[0]
    rho_list=[rho_update]
    convergence='False'
    while convergence=='False':
        i+=1
        print('iteration '+str(i)+'...')
        step_list+=[i]
        previous_LL=Init_LL[-1]
        t1t2_list=OPTIMIZE_t1t2(overlapping_arrays, rho_update, size_lims)
        (rho_update, neg_LL_update) = OPTIMIZE_rho(t1t2_list,pair_list,size_lims,non_overlapping_arrays,ancestor_dicts)
        Init_LL+=[neg_LL_update[0]]
        rho_list+=[rho_update]
        if neg_LL_update[0]>previous_LL:
            convergence='True'
    final_dist=dict(zip(pair_list,t1t2_list))
    return(final_dist,rho_list[-1:])

######### Phylogeny from CRISPR_distance ###########

def get_weights_for_all_pairs(df):
    '''
    Given a non-reversible distance matrix as pandas DataFrame, compute the weights for all pairs (x, y). Outputs the pair with maximum weight
    '''
    wt=pd.DataFrame(np.nan,columns=df.columns,index=df.index)

    for tu in list(itertools.permutations(df.columns,2)):
        x=tu[0]
        y=tu[1]
        d_xy=df.loc[x,y]
        d_yx=df.loc[y,x]
        d_xz=sum(df.loc[x,:].drop([x,y]).fillna(0))
        d_yz=sum(df.loc[y,:].drop([x,y]).fillna(0))
        wt.loc[x,y]=((2-len(df.columns))*(d_xy+d_yx))+d_xz+d_yz
    
    max_x=wt.max().idxmax()
    max_y=wt.loc[wt.max().idxmax()].idxmax()
    max_pair=(max_x,max_y)
    return(max_pair)

def get_new_node_distance(df,max_pair):
    ''' Given a non-reversible distance matrix as pandas DataFrame and the pair of closest neighbors, create new node and calculate distance to new node. Outputs the non-reversible distance matric with new node. Note that the initial pair is not dropped from the df  '''
    df2=df.copy()
    t=df.index.max()+1
    x=max_pair[0]
    y=max_pair[1]
    df2.loc[t] = ((df.loc[x,:]-df.loc[x,y]).fillna(0)+(df.loc[y,:]-df.loc[y,x]).fillna(0))/2
    df2[t] = (df.loc[:,x].fillna(0)+df.loc[:,y].fillna(0))/2
    df2.loc[x,t]=df.loc[x,y]
    df2.loc[y,t]=df.loc[y,x]
    return(df2)

def phylogeny_from_CRISPR(arrays,final_dist):
    '''Given a list of CRISPR arrays and a dictionnary of optimized distances, returns a phylogenetic tree (Tree class) of the CRISPR arrays and a dictionary of labels'''
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

    # init terminal clades 
    clades = dict(zip(df.index,(BaseTree.Clade(None, str(name)) for name in df.index)))
    inner_count = len(df)-1

    while len(df.columns)>1:
        inner_count+=1

        # find minimum distance pair
        max_pair=get_weights_for_all_pairs(df)
        min_i=max_pair[0]
        min_j=max_pair[1]

        # calculate nodeDist
        df=get_new_node_distance(df,max_pair)

        # create clade 
        clade1 = clades[min_i] 
        clade2 = clades[min_j] 
        inner_clade = BaseTree.Clade(None, str(inner_count)) 
        inner_clade.clades.append(clade1) 
        inner_clade.clades.append(clade2) 
        # assign branch length 
        clade1.branch_length = float(df.loc[min_i, min_j])
        clade2.branch_length = float(df.loc[min_j, min_i])
        # update node dict 
        clades[inner_count] = inner_clade 
        del clades[min_i]
        del clades[min_j]
        if len(clades)==1:
            inner_clade.branch_length = None

        # rebuild distance matrix
        df=df.drop([min_i,min_j],axis=1).drop([min_i,min_j],axis=0)

    TREE=Tree(inner_clade,rooted=True)
    return(TREE,dic)


def get_distance_matrix_from_phylogeny(Tree):
    '''
    Given a phylogenetic tree Tree, returns a reversible euclidian distance matrix (as a pandas dataframe)
    '''
    taxa=sorted(map(int,[cl.name for cl in Tree.get_terminals()]))
    matrix=pd.DataFrame([],columns=taxa,index=taxa)
    for combi in itertools.combinations(range(len(taxa)),2):
        tgt1=str(combi[0])
        tgt2=str(combi[1])
        dist=Tree.distance({"name": tgt1}, {"name": tgt2})
        matrix.loc[combi]=dist
        matrix.loc[combi[1],combi[0]]=dist
    return(matrix)

############ PARALLELIZED VERSION ################

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
    print(f'number of processeors: {nproc}')
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
    print('\nGenerating ancestors...')
    argList = [tuple([array,size_lims]) for array in overlapping_arrays]
    with multiprocessing.Pool(processes=nproc) as pool:
        ancestor_dicts = pool.map(get_ancestor_dict,argList)
    print('Complete.\n')

    # first iteration: get initial t1s,t2s, and rho
    print('Initial estimates...')
    with multiprocessing.Pool(processes=nproc) as pool:
        # create argList for optimization of t1 and t2
        argList = [(pair,rho_init,size_lims,ancestor_dicts[index])
                   for index,pair in enumerate(pairList)]
        print('Estimating initial t1 and t2 lists...')
        t1t2_list = pool.map(OPTIMIZE_pair_t1t2,argList)
    print('Complete. Updating rho...')
    rho_update, neg_LL_update = OPTIMIZE_rho(t1t2_list,pairList,size_lims,non_overlapping_arrays,ancestor_dicts)
    print('Complete. Starting main interations.\n')
    
    Init_LL=[np.inf]
    i=0
    j=0
    max_iters_after_converged = 5
    rho_list=[rho_update[0]]
    convergence=False
    
    # optimization of rho, and (t1, t2) of each pair
    while ((convergence==False) or (j <= max_iters_after_converged)):
        i+=1
        print('Iteration '+str(i)+'...')
        previous_LL=Init_LL[-1]
        # optimize t1 and t2
        with multiprocessing.Pool(processes=nproc) as pool:
            argList = [(pair,rho_update,size_lims,ancestor_dicts[index])
                       for index,pair in enumerate(pairList)]
            t1t2_list = pool.map(OPTIMIZE_pair_t1t2,argList)
        # optimize rho and return neg log likelihood
        rho_update, neg_LL_update = OPTIMIZE_rho(t1t2_list,pairList,size_lims,non_overlapping_arrays,ancestor_dicts)
        # check for convergence
        Init_LL+=[neg_LL_update[0]]
        rho_list+=[rho_update[0]]
        if neg_LL_update[0]>previous_LL:
            convergence=True
        if convergence == True: # iterating beyond convergence
            j+=1
    print('Converged.')
    final_dist=dict(zip(pairList,t1t2_list))
    LOG = pd.DataFrame([Init_LL,rho_list],index=['LogLikelihood','Rho']).T

    # Get phylogeny and distance matrix
    Tree,labels=phylogeny_from_CRISPR(arrays,final_dist)
    dist=get_distance_matrix_from_phylogeny(Tree)
    print(f'Final rho: {rho_list[-1]}')
    return(Tree,dist,labels,LOG,final_dist)