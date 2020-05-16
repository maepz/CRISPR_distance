#!/usr/bin/env python
# coding: utf-8

# # Here are the functions needed for CRISPR_distance

# ## Get the ansestors
def get_limits_ancestor_sizes(arrays,p_low=0.005,p_high=0.995):
    '''Given a set of arrays and rho, returns l1 and l2 the min and max ancestor lengths'''

    from scipy.stats import poisson
    rho=sum([len(ar) for ar in arrays])/len([len(ar) for ar in arrays])
    return(tuple(poisson.ppf([p_low,p_high],rho).astype(int)))

def nCr(n,r):
        import math
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
        from CRISPR_functions import nCr
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
        d1=len(PAIR.d1_spacers)
        d2=len(PAIR.d2_spacers)


        m1=c1+i
        d1=j+u+c-c1
        j1=d1-i
        m2=c2+j
        d2=i+u+c-c2
        j2=d2-j

        return(m1,d1,j1,m2,d2,j2)
        
    def get_arrays(self,c_spacers,d1_spacers,d2_spacers): # to get the list of all possible arrays not usefull       
        array_list2=[]
        def merge_lists(lst1,lst2):
            import itertools
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
        from CRISPR_functions import categorize_spacers_for_ordered_model
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
        from CRISPR_functions import combi
       
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
            ancestor_weights = [combi.array_counts/ancestor_counts for combi in v]
            
            combi_output.update(dict(zip(combi_strings,ancestor_weights)))

        return(combi_output)


# ## Math equations

### Length model 

def rho(lbda,mu):
    import numpy as np
    return(lbda/mu)

def prob_n_given_ro(n,rho):
    import numpy as np
    return((np.e**-rho)*(rho**n/np.math.factorial(n))) 

### Transitions prob for the ordered model

def M(m,t,mu):
    '''probability of preserving m spacers in time t : M(m|t,mu)'''
    import numpy as np
    M=np.exp(-m*t*mu)
    return(M)

def D(d,t,mu):
    '''probability of loosing d spacers in time t : D(d|t,mu)'''
    import numpy as np
    D=np.power((1-np.exp(-t*mu)),d)
    return(D)

def I(j,t,lbda,mu):
    '''probability of inserting d spacers in time t : D(d|t,mu)'''
    import numpy as np
    rho=lbda/mu
    a=rho*(1-np.exp(-t*mu))
    I=(np.exp(-a)*np.power(a,j))/np.math.factorial(j)
    return(I)


def L(rho,t,PAIR,size_lims):
    from CRISPR_functions import prob_n_given_ro, combi, CRISPR_pair, M, I, D
    t1=t[0]
    t2=t[1]
    l1=size_lims[0]
    l2=size_lims[1]
    lbda=rho/(rho+1)
    mu=1/(rho+1)
    
    ancestors_dic= PAIR.get_combi(size_lims)
    Likelihood=0
    
    for s0,ws in ancestors_dic.items():
        combi_s0 = combi(list(map(int,s0.split('-'))))
        n=combi_s0.n

        m1,d1,j1,m2,d2,j2 = combi_s0.compare_ancestor_to_arrays(PAIR)

        Qa = prob_n_given_ro(n,rho)*ws
        T1 = M(m1,t1,mu)*D(d1,t1,mu)*I(j1,t1,lbda,mu)
        T2 = M(m2,t2,mu)*D(d2,t2,mu)*I(j2,t2,lbda,mu)
        
        Likelihood += Qa*T1*T2
    
    return(Likelihood)

def neg_LL_floating_t(x,rho,PAIR,size_lims):
    import numpy as np
    return(-np.log(L(rho=rho,t=x,PAIR=PAIR,size_lims=size_lims)))
    
def neg_LL_floating_rho(x,t1t2_list,pair_list,size_lims,non_overlapping_arrays):
    import numpy as np
    neg_LL_overlapping=sum([-np.log(L(rho=x,t=t1t2_list[i],PAIR=pair_list[i],size_lims=size_lims)) for i in range(len(pair_list))])
    neg_LL_non_overlapping=sum([-np.log(prob_n_given_ro(len(pair[0]),x)*prob_n_given_ro(len(pair[1]),x)) for pair in non_overlapping_arrays])
    return(neg_LL_overlapping+neg_LL_non_overlapping)


def OPTIMIZE_rho(t1t2_list,pair_list,size_lims,non_overlapping_arrays):
    from scipy.optimize import minimize
    return(float(minimize(neg_LL_floating_rho,[2],method='Powell',args=(t1t2_list,pair_list,size_lims,non_overlapping_arrays)).x))

def OPTIMIZE_t1t2(overlapping_arrays, rho, size_lims):
    '''Provided a set of overlapping arrays and rho, find their best respective divergence times t1,t2 from ancestor to arrays. The output is an array of [t1,t2] of length len(overlapping_arrays)'''
    from scipy.optimize import minimize
    
    t1t2_list=[]

    for pair in overlapping_arrays:
        PAIR=CRISPR_pair(pair[0],pair[1])
        t1t2_list+=[tuple(minimize(neg_LL_floating_t,[1,1],method='Powell',args=(rho,PAIR,size_lims)).x)]
    
    return(t1t2_list)


