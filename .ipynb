{
 "cells": [
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "# Introduction\n",
    "This script should calculate a pairwise distance matrix based on the CRISPR arrays.\n",
    "It is an adaptation of the method in Kupczok and Bollack 2013*\n",
    "\n",
    "*Kupczok, A., and Bollback, J.P. 2013. Probabilistic models for CRISPR spacer content evolution. BMC Evolutionary Biology 13(1): 54. doi:10.1186/1471-2148-13-54.\n"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## From the paper\n",
    "1. Estimate a starting ρ from the length model by maximum likelihood. \n",
    "2. For each pair of spacers with overlap, generate the\n",
    "    possible ancestors: Ancestral arrays can be arbitrarily\n",
    "    large, but the probability of observing a certain\n",
    "    length is given by p(n). For practical reasons we do\n",
    "    not consider ancestors whose length is outside the\n",
    "    central 99% of the stationary distribution given by ρ\n",
    "    estimated in step 1, since they would have a\n",
    "    negligible contribution to the likelihood. In detail, the\n",
    "    length l1 where the cumulative distribution exceeds\n",
    "    0.005 is the minimum ancestor length and the length\n",
    "    l2 where the cumulative distribution exceeds 0.995 is\n",
    "    the maximum ancestor length. Then the possible\n",
    "    ancestor lengths n are between l1 and l2: l1 ≤ n ≤ l2.\n",
    "3. (a) For all pairs with overlap, estimate the times\n",
    "    with fixed ρ. It is possible to iterate through\n",
    "    the pairs and estimate their times\n",
    "    independently of the other pairs. The\n",
    "    estimation of both times is iterated\n",
    "    alternatingly until the likelihood has\n",
    "    converged.\n",
    "    (b) Estimate ρ with fixed times using L(ρ|t, S).\n",
    "    (c) Check if the log-likelihood of the estimated\n",
    "    parameters has converged, then return the\n",
    "    estimated parameters, else repeat step (a)\n",
    "    with the new parameters."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Pseudo code\n",
    "\n",
    "functions:\n",
    "Amongst the spacers, find each pair of spacers with overlap\n",
    "Generate the possible ancestors of a pair\n",
    "Estimate time divergence given fixed p \n",
    "Estimate p given fixed divergence time\n"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "# Find pairs with overlap"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 161,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "1"
      ]
     },
     "execution_count": 161,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "ancestor=[0,1,2,3,4,5,6,7,8,9,10,11,12,13]\n",
    "arr=[[9,2,3,4,5],[0,1,2,3,7,8],[1,10,11,12,13]]\n",
    "# Deconstruct array into fragments\n",
    "array1=arr[0]\n",
    "array2=arr[1]    \n",
    "\n",
    "def is_overlapping(array1, array2):\n",
    "    for i in array1:\n",
    "        if i in array2:\n",
    "            answer=1\n",
    "            break\n",
    "        else:\n",
    "            answer=0\n",
    "            continue\n",
    "    return(answer)\n",
    "\n",
    "is_overlapping(array1, array2)    "
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "# "
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Ordered model:\n",
    "Given a pair of arrays s1 and s2, find the\n",
    "first shared spacer. The ancestor must contain this spacer\n",
    "and all subsequent spacers from both arrays, these are c\n",
    "spacers in total. There are d1 and d2 spacers before the\n",
    "first shared spacer in s1 and s2, respectively. With these\n",
    "new definitions of c, d1 and d2, the method from the\n",
    "unordered model is applied."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 419,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "s1 [9, 2, 3, 4, 5]\n",
      "s2 [0, 1, 2, 3, 7, 8]\n",
      "c_spacers [2, 3, 4, 5, 7, 8]\n",
      "d1_spacers [9]\n",
      "d2_spacers [0, 1]\n"
     ]
    }
   ],
   "source": [
    "'''Given a pair of arrays s1 and s2, find the\n",
    "first shared spacer. The ancestor must contain this spacer\n",
    "and all subsequent spacers from both arrays, these are c\n",
    "spacers in total. There are d1 and d2 spacers before the\n",
    "first shared spacer in s1 and s2, respectively. With these\n",
    "new definitions of c, d1 and d2, the method from the\n",
    "unordered model is applied.'''\n",
    "\n",
    "def categorize_spacers_for_ordered_model(s1,s2):\n",
    "    c_spacers=[]\n",
    "    for i in range(len(s1)):\n",
    "        sp=s1[i]\n",
    "        if sp in s2:\n",
    "            j=array2.index(sp)\n",
    "            c_spacers=s2[j:]+s1[i:]\n",
    "            d1_spacers=s1[:i]\n",
    "            d2_spacers=s2[:j]\n",
    "            break\n",
    "    return(list(set(c_spacers)),d1_spacers,d2_spacers)\n",
    "\n",
    "arr=[[9,2,3,4,5],[0,1,2,3,7,8],[1,10,11,12,13]]\n",
    "s1=arr[0]\n",
    "s2=arr[1]\n",
    "print('s1',s1)\n",
    "print('s2',s2)\n",
    "#first spaer of shorter array\n",
    "\n",
    "c_spacers, d1_spacers, d2_spacers = categorize_spacers_for_ordered_model(s1,s2)           \n",
    "print('c_spacers',c_spacers)\n",
    "print('d1_spacers',d1_spacers)\n",
    "print('d2_spacers',d2_spacers)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "Then all n between min(c, l1) and l2 are generated.\n",
    "When length n is generated, enumerate all i, j, u\n",
    "such that c + i + j + u = n, i ≤ d1 and j ≤ d2. Then for\n",
    "ancestor a, there are c common spacers, i only occur in s1, j\n",
    "only occur in s2 and u are unobserved (they are lost in both\n",
    "lineages)."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "# Independent loss models : Length model\n",
    "ro = $$\\rho = \\frac{\\lambda}{\\mu}$$\n",
    "with $\\lambda$ = spacer insertion rate, $\\mu$ = spacer_deletion_rate\n",
    "\n",
    "prob_n_given_ro = $$p(n|\\rho) = e^{-\\rho} \\frac{\\rho^n}{n!}$$"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 412,
   "metadata": {},
   "outputs": [],
   "source": [
    "import numpy as np\n",
    "\n",
    "lbda=1\n",
    "mu=2\n",
    "n=2\n",
    "\n",
    "def rho(lbda,mu):\n",
    "    return(lbda/mu)\n",
    "def prob_n_given_ro(n,rho):\n",
    "    return((np.e**-rho)*(rho**n/np.math.factorial(n))) "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 360,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "36.0\n",
      "36\n"
     ]
    }
   ],
   "source": [
    "#### define combi class ####\n",
    "class combi:\n",
    "    \n",
    "    \"\"\"The combi object is a list of integers [c,i,j,u] that represent all putative ancestors given the following parameters:\n",
    "    c = number of c_spacers; spacers necessarily present in the ancestor of two CRISPR arrays,\n",
    "    i = number of putative ancestral spacers amongst the spacers only present in array1,\n",
    "    j = number of putative ancestral spacers amongt these only present in array 2,\n",
    "    u = number of putative spacers present in ancestor but lost in both arrays\n",
    "    \n",
    "    the lengths of the putative ancestral arrays n = sum([c,i,j,u])\"\"\"\n",
    "    def nCr(n,r):\n",
    "        import math\n",
    "        f = math.factorial\n",
    "        return f(n) / f(r) / f(n-r)\n",
    "    \n",
    "    def __init__(self,list):\n",
    "        self.c= list[0]\n",
    "        self.i= list[1]\n",
    "        self.j= list[2]\n",
    "        self.u= list[3]\n",
    "        self.n= sum(list)\n",
    "        possible_positions_of_u=nCr(self.n,self.u)\n",
    "        possible_combinations_of_d1_and_d2=nCr(self.i+self.j,self.i)\n",
    "        \n",
    "        self.array_counts=possible_positions_of_u*possible_combinations_of_d1_and_d2 # to get the number of putative ancestors\n",
    "        \n",
    "    def get_arrays(self,c_spacers,d1_spacers,d2_spacers): # to get the list of all possible arrays        \n",
    "        array_list2=[]\n",
    "        def merge_lists(lst1,lst2):\n",
    "            import itertools\n",
    "            array_list=[]\n",
    "            for locations in itertools.combinations(range(len(lst1) + len(lst2)), len(lst2)):\n",
    "                result = lst1[:]\n",
    "                for location, element in zip(locations, lst2):\n",
    "                    result.insert(location, element)\n",
    "                new_list=result\n",
    "                array_list+=[new_list]\n",
    "            return(array_list)\n",
    "        \n",
    "        d1_good=d1_spacers[len(d1_spacers)-self.i:]\n",
    "        d2_good=d2_spacers[len(d2_spacers)-self.j:]\n",
    "        u_good=self.u*['u']\n",
    "        c_spacers=c_spacers\n",
    "        if self.u>0:\n",
    "            for array in merge_lists(d1_good,d2_good):\n",
    "#                 print (u_good,array+c_spacers)\n",
    "                array_list2+=merge_lists(u_good,array+c_spacers)\n",
    "        else:\n",
    "            array_list2=[array+c_spacers for array in merge_lists(d1_good,d2_good)]\n",
    "        return(array_list2)\n",
    "        \n",
    "# print(combi([6, 0, 2, 1]).c)\n",
    "print(combi([6, 1, 0, 2]).array_counts)\n",
    "print(len(combi([6, 1, 0, 2]).get_arrays(c_spacers,d1_spacers,d2_spacers)))\n",
    "# combi([6, 1, 0, 2]).get_arrays(c_spacers,d1_spacers,d2_spacers)\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 2,
   "metadata": {},
   "outputs": [
    {
     "ename": "ModuleNotFoundError",
     "evalue": "No module named 'CRISPR_distance'",
     "output_type": "error",
     "traceback": [
      "\u001b[0;31m---------------------------------------------------------------------------\u001b[0m",
      "\u001b[0;31mModuleNotFoundError\u001b[0m                       Traceback (most recent call last)",
      "\u001b[0;32m<ipython-input-2-b9d30c59c9e6>\u001b[0m in \u001b[0;36m<module>\u001b[0;34m\u001b[0m\n\u001b[1;32m     24\u001b[0m \u001b[0ms2\u001b[0m\u001b[0;34m=\u001b[0m\u001b[0marr\u001b[0m\u001b[0;34m[\u001b[0m\u001b[0;36m1\u001b[0m\u001b[0;34m]\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m     25\u001b[0m \u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0;32m---> 26\u001b[0;31m \u001b[0mpair\u001b[0m\u001b[0;34m=\u001b[0m\u001b[0mCRISPR_pair\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0ms1\u001b[0m\u001b[0;34m,\u001b[0m\u001b[0ms2\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0m\u001b[1;32m     27\u001b[0m \u001b[0mpair\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0mc_spacers\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n",
      "\u001b[0;32m<ipython-input-2-b9d30c59c9e6>\u001b[0m in \u001b[0;36m__init__\u001b[0;34m(self, s1, s2)\u001b[0m\n\u001b[1;32m     15\u001b[0m \u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m     16\u001b[0m     \u001b[0;32mdef\u001b[0m \u001b[0m__init__\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mself\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0ms1\u001b[0m\u001b[0;34m,\u001b[0m\u001b[0ms2\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m:\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0;32m---> 17\u001b[0;31m         \u001b[0;32mfrom\u001b[0m \u001b[0mCRISPR_distance\u001b[0m \u001b[0;32mimport\u001b[0m \u001b[0mcategorize_spacers_for_ordered_model\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0m\u001b[1;32m     18\u001b[0m         \u001b[0mself\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0ms1\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0ms1\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m     19\u001b[0m         \u001b[0mself\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0ms2\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0ms2\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n",
      "\u001b[0;31mModuleNotFoundError\u001b[0m: No module named 'CRISPR_distance'"
     ]
    }
   ],
   "source": [
    "class CRISPR_pair:\n",
    "    \n",
    "    def categorize_spacers_for_ordered_model(s1, s2):\n",
    "        c_spacers=[]\n",
    "        \n",
    "        for i in range(len(self.s1)):\n",
    "            sp=s1[i]\n",
    "            if sp in s2:\n",
    "                j=array2.index(sp)\n",
    "                c_spacers=s2[j:]+s1[i:]\n",
    "                d1_spacers=s1[:i]\n",
    "                d2_spacers=s2[:j]\n",
    "                break\n",
    "        return(list(set(c_spacers)),d1_spacers,d2_spacers)\n",
    "    \n",
    "    def __init__(self, s1,s2):\n",
    "        from CRISPR_distance import categorize_spacers_for_ordered_model\n",
    "        self.s1 = s1\n",
    "        self.s2 = s2\n",
    "        self.c_spacers, self.d1_spacers, self.d2_spacers = categorize_spacers_for_ordered_model(self.s1,self.s2)\n",
    "    \n",
    "arr=[[9,2,3,4,5],[0,1,2,3,7,8],[1,10,11,12,13]]\n",
    "s1=arr[0]\n",
    "s2=arr[1]\n",
    "\n",
    "pair=CRISPR_pair(s1,s2)\n",
    "pair.c_spacers"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 397,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "9\n",
      "1.0\n",
      "9.0\n",
      "47.0\n",
      "186.0\n"
     ]
    },
    {
     "data": {
      "text/plain": [
       "{6: {(6, 0, 0, 0): {'w': 1.0, 'ws': 1.0}},\n",
       " 7: {(6, 0, 0, 1): {'w': 7.0, 'ws': 0.7777777777777778},\n",
       "  (6, 0, 1, 0): {'w': 1.0, 'ws': 0.1111111111111111},\n",
       "  (6, 1, 0, 0): {'w': 1.0, 'ws': 0.1111111111111111}},\n",
       " 8: {(6, 0, 0, 2): {'w': 28.0, 'ws': 0.5957446808510638},\n",
       "  (6, 0, 1, 1): {'w': 8.0, 'ws': 0.1702127659574468},\n",
       "  (6, 0, 2, 0): {'w': 1.0, 'ws': 0.02127659574468085},\n",
       "  (6, 1, 0, 1): {'w': 8.0, 'ws': 0.1702127659574468},\n",
       "  (6, 1, 1, 0): {'w': 2.0, 'ws': 0.0425531914893617}},\n",
       " 9: {(6, 0, 0, 3): {'w': 84.0, 'ws': 0.45161290322580644},\n",
       "  (6, 0, 1, 2): {'w': 36.0, 'ws': 0.1935483870967742},\n",
       "  (6, 0, 2, 1): {'w': 9.0, 'ws': 0.04838709677419355},\n",
       "  (6, 1, 0, 2): {'w': 36.0, 'ws': 0.1935483870967742},\n",
       "  (6, 1, 1, 1): {'w': 18.0, 'ws': 0.0967741935483871},\n",
       "  (6, 1, 2, 0): {'w': 3.0, 'ws': 0.016129032258064516}}}"
      ]
     },
     "execution_count": 397,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "#### get w and ws ####\n",
    "\n",
    "c=len(c_spacers)\n",
    "# print(len(d1_spacers))\n",
    "n=len(c_spacers+d1_spacers+d2_spacers)\n",
    "len(d1_spacers)\n",
    "len(d2_spacers)\n",
    "# n=len(ancestor)\n",
    "# n=[n for n in range(min(len(c_spacers),len(s1)),len(ancestor))]\n",
    "print(n)\n",
    "''' Get the possible combinations of spacers categories for each ancestor length. n:[[c,i,j,u]] with n length of ancestral array, c number of spacers in common (spacers necessarily present in ancerstor), i number of ancestral spacers amongst the spacers only present in array1, j number of ancestral spacers amongt these only present in array 2'''\n",
    "spacers_combi={}\n",
    "for n in range(c,n+1):\n",
    "#     print('n',n)\n",
    "    i=0\n",
    "    spacers_combi[n]={}\n",
    "    for i in range(len(d1_spacers)+1):\n",
    "        for j in range(len(d2_spacers)+1):\n",
    "            if c+i+j<n+1:\n",
    "#                 print('j',j)\n",
    "\n",
    "                u=n-(c+i+j)\n",
    "#                 print(c,i,j,u)\n",
    "                spacers_combi[n][c,i,j,u]={'w':combi([c,i,j,u]).array_counts}\n",
    "    \n",
    "                    \n",
    "for k,v in spacers_combi.items():\n",
    "    ktot=0\n",
    "    for val in v.values():\n",
    "        ktot+=val['w']\n",
    "    print(ktot)\n",
    "    for c in v.keys():\n",
    "        spacers_combi[k][c].update({'ws':spacers_combi[k][c]['w']/ktot})\n",
    "spacers_combi"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 423,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "2.0"
      ]
     },
     "execution_count": 423,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "nCr(len(d1_spacers),0)*nCr(len(d2_spacers),1)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 356,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "[2, 3, 4, 5, 7, 8]\n",
      "[9]\n",
      "[0, 1]\n",
      "combi (6, 0, 0, 3)\n",
      "84.0\n",
      "combi (6, 0, 1, 2)\n",
      "36.0\n",
      "combi (6, 0, 2, 1)\n",
      "9.0\n",
      "combi (6, 1, 0, 2)\n",
      "36.0\n",
      "combi (6, 1, 1, 1)\n",
      "18.0\n",
      "combi (6, 1, 2, 0)\n",
      "3.0\n"
     ]
    }
   ],
   "source": [
    "#### get array counts (now within combi class) ####\n",
    "print(c_spacers)\n",
    "print(d1_spacers)\n",
    "print(d2_spacers)\n",
    "n=9\n",
    "for combi in spacers_combi[n].keys():\n",
    "    print ('combi',combi)\n",
    "    c=combi[0]\n",
    "    i=combi[1]\n",
    "    j=combi[2]\n",
    "    u=combi[3]\n",
    "    possible_positions_of_u=nCr(n,u)\n",
    "    possible_combinations_of_d1_and_d2=nCr(i+j,i)\n",
    "    arrays_for_this_combi=possible_positions_of_u*possible_combinations_of_d1_and_d2\n",
    "    print(arrays_for_this_combi)\n",
    "#     d1_good=d1_spacers[len(d1_spacers)-i:]\n",
    "#     d2_good=d2_spacers[len(d2_spacers)-j:]\n",
    "#     u_good=u*'u'\n",
    "#     print(arrays_for_this_combi)\n",
    "#     print(d1_good)\n",
    "#     print(d2_good)\n",
    "    "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 409,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "[[0, 1, 9, 2, 3, 4, 5, 7, 8],\n",
       " [0, 9, 1, 2, 3, 4, 5, 7, 8],\n",
       " [9, 0, 1, 2, 3, 4, 5, 7, 8]]"
      ]
     },
     "execution_count": 409,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "#### get ancestor arrays (now withon combi class) ####\n",
    "\n",
    "def get_ancestor_arrays(c,i,j,u,d1_spacers,d2_spacers,c_spacers):\n",
    "    array_list=[]\n",
    "    def merge_lists(lst1,lst2):\n",
    "        import itertools\n",
    "        array_list=[]\n",
    "\n",
    "        for locations in itertools.combinations(range(len(lst1) + len(lst2)), len(lst2)):\n",
    "            result = lst1[:]\n",
    "            for location, element in zip(locations, lst2):\n",
    "                result.insert(location, element)\n",
    "            new_list=result\n",
    "            array_list+=[new_list]\n",
    "        return(array_list)\n",
    "    \n",
    "    d1_good=d1_spacers[len(d1_spacers)-i:]\n",
    "    d2_good=d2_spacers[len(d2_spacers)-j:]\n",
    "    u_good=u*'u'\n",
    "    c_spacers=c_spacers\n",
    "    if u>0:\n",
    "        for array in merge_lists(d1_good,d2_good):\n",
    "            print ([u_good],array+c_spacers)\n",
    "            array_list+=merge_lists(u_good,array+c_spacers)\n",
    "    else:\n",
    "        array_list=[array+c_spacers for array in merge_lists(d1_good,d2_good)]\n",
    "    return(array_list)\n",
    "    \n",
    "get_ancestor_arrays(c,i,j,u,d1_good,d2_good,c_spacers)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 2
}
