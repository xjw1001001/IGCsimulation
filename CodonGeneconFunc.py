# Xiang decided to seperate functions and the recording class so that he can easily vectorize his code
# xji3@ncsu.edu

# Uses Alex Griffing's JsonCTMCTree package for likelihood and gradient calculation
# Re-written of my previous CondonGeneconv class

from operator import mul
from itertools import product
from functools import partial
from copy import deepcopy
import pickle
import os, sys

import numpy as np
import networkx as nx
import scipy
import scipy.optimize
import scipy.sparse
import scipy.sparse.linalg

from Bio import Phylo
from Bio import SeqIO

import cProfile
import jsonctmctree.ll, jsonctmctree.interface

##def get_HKYBasicRate(na, nb, pi, kappa):
##    if isTransition(na, nb):
##        return pi['ACGT'.index(nb)] * kappa
##    else:
##        return pi['ACGT'.index(nb)]

def get_HKYGeneconvRate(pair_from, pair_to, Qbasic, tau):
    na, nb = pair_from
    nc, nd = pair_to
    if (na != nc and nb!= nd) or (na == nc and nb == nd):
        return 0.0
    if na ==nc and nb != nd:
        Qb = Qbasic['ACGT'.index(nb), 'ACGT'.index(nd)]
        if na == nd:
            return Qb + tau
        else:
            return Qb
    if nb == nd and na != nc:
        Qb = Qbasic['ACGT'.index(na), 'ACGT'.index(nc)]
        if nb == nc:
            return Qb + tau
        else:
            return Qb
    print ('Warning: Check get_HKYGeneconvRate Func. You should not see this.')

def get_MG94BasicRate(ca, cb, pi, kappa, omega, codon_table):
    dif = [ii for ii in range(3) if ca[ii] != cb[ii]]
    ndiff = len(dif)
    if ndiff > 1:
        return 0
    elif ndiff == 0:
        print 'Please check your codon tables and make sure no redundancy'
        print ca, cb
        return 0
    else:
        na = ca[dif[0]]
        nb = cb[dif[0]]
        QbasicRate = pi['ACGT'.index(nb)]

        if isTransition(na, nb):
            QbasicRate *= kappa

        if isNonsynonymous(ca, cb, codon_table):
            QbasicRate *= omega

        return QbasicRate

def isTransition(na, nb):
    return (set([na, nb]) == set(['A', 'G']) or set([na, nb]) == set(['C', 'T']))

def isNonsynonymous(ca, cb, codon_table):
    return (codon_table[ca] != codon_table[cb])

#vec_get_MG94BasicRate = np.vectorize(get_MG94BasicRate, doc='Vectorized `get_MG94BasicRate`', excluded = ['pi', 'kappa', 'omega', 'codon_table'])

def get_MG94GeneconvRate(pair_from, pair_to, Qbasic, tau, codon_to_state):
    # pair_from = a string of length 6
    ca, cb = pair_from[:3], pair_from[3:]
    cc, cd = pair_to[:3], pair_to[3:]
    row = (codon_to_state[ca], codon_to_state[cb])
    col = (codon_to_state[cc], codon_to_state[cd])
    if ca != cc and cb != cd:
        return None
    if ca == cc and cb == cd: # do not deal with diagonal entries here
        return None

    if ca == cc:
        # cb != cd
        BasicRate = Qbasic[codon_to_state[cb], codon_to_state[cd]]
        if cd == ca:
            if isNonsynonymous(cb, cd, codon_table):
                additional_source = tau * omega
            else:
                additional_source = tau
            return [row, col, BasicRate + additional_source]
        else:
            return [row, col, BasicRate]
    else: # cb == cd
        # ca != cc
        BasicRate = Qbasic[codon_to_state[ca], codon_to_state[cc]]
        if cc == cb:
            if isNonsynonymous(ca, cc, codon_table):
                additional_source = tau * omega
            else:
                additional_source = tau
            return [row, col, BasicRate + additional_source]
        else:
            return [row, col, BasicRate]
    #print ('You should not see this')

## vec_get_MG94GeneconvRate = np.vectorize(get_MG94GeneconvRate, doc = 'Vectorized `get_MG94GeneconvRate`', excluded = ['pi', 'kappa', 'omega', 'codon_table', 'tau', 'codon_to_state'])

def get_x_clock_guess(edge_to_blen):
    # TODO: modify this to work on general tree topology
    leaf_branch = [edge for edge in edge_to_blen.keys() if edge[0][0] == 'N' and str.isdigit(edge[0][1:]) and not str.isdigit(edge[1][1:])]
    out_group_branch = [edge for edge in leaf_branch if edge[0] == 'N0' and not str.isdigit(edge[1][1:])] [0]
    internal_branch = [x for x in edge_to_blen.keys() if not x in leaf_branch]


    leaf_branch.sort(key = lambda node: int(node[0][1:]))  # sort the list by the first node number in increasing order
    internal_branch.sort(key = lambda node: int(node[0][1:]))  # sort the list by the first node number in increasing order

    Lr_reverse = []
    for i in range(len(internal_branch) - 1, 0, -1):
        Lr_reverse.append(edge_to_blen[leaf_branch[i]]/(edge_to_blen[leaf_branch[i]] + edge_to_blen[internal_branch[i]]))
    r0 = 2 / (edge_to_blen[out_group_branch] * (1 - Lr_reverse[-1]) / edge_to_blen[internal_branch[0]] + 1)
    Lr_reverse.append(r0)
    Lr_reverse.append(edge_to_blen[out_group_branch] / (2 - r0))
    return Lr_reverse

if __name__ == '__main__':
    print 
                                        
    
