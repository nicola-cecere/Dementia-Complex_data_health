#! /usr/bin/env python

"""
# -----------------------------------------------------------------------
#
# localization.py
#
# by Joerg Menche
# ADAPTED VERSION (2023-11-06)
# Last Modified: 2014-12-06
#
# This code determines the network-based location in a module for a given
# set of nodes on a network as described in 
# 
# Uncovering Disease-Disease Relationships Through The Human
# Interactome
#
# by Joerg Menche, Amitabh Sharma, Maksim Kitsak, Susan Dina
#    Ghiassian, Marc Vidal, Joseph Loscalzo & Albert-Laszlo Barabasi
# 
# -----------------------------------------------------------------------
"""


import networkx as nx
import random 
import numpy as np
from scipy import stats


# =================================================================================
def get_lcc(G,seed_nodes):
    
    g = nx.subgraph(G, seed_nodes)
    
    lcc = max(nx.connected_components(g), key=len)

    return lcc


# =============================================================================
def get_lcc_size(G,seed_nodes):

    # getting subgraph that only consists of the black_nodes
    g = nx.subgraph(G,seed_nodes)

    lcc = max(nx.connected_components(g), key=len)
        
    return len(lcc)


# =============================================================================
def group_nodes_by_degree(G):

    degree_dict = {}

    for node in G.nodes():
        degree = G.degree(node)

        if degree not in degree_dict:
            degree_dict[degree] = []

        degree_dict[degree].append(node)

    return degree_dict


# =============================================================================
def sample_preserving_degrees(G,S,bucket):
    
    sampled_nodes = set()

    for node in S:
        degree = G.degree(node)

        available_nodes = bucket[degree]

        flag = True
        while(flag):
            chosen_node = random.choice(available_nodes)
            flag = chosen_node in sampled_nodes

        sampled_nodes.add(chosen_node)
    return sampled_nodes


# =============================================================================
def get_random_comparison(G,genes,sims,degree_preserving=False, hist=None):

    # getting all genes in the network  
    all_genes = sorted(G.nodes())

    gene_set = set(genes) & set(G.nodes())

    number_of_seed_genes = len(gene_set & set(all_genes))
    
    l_list  = []

    if degree_preserving:
        bucket = group_nodes_by_degree(G)
        

    # simulations with randomly distributed seed nodes
    for i in range(sims):
        if degree_preserving:
            rand_seeds = sample_preserving_degrees(G, gene_set, bucket)
        else:
            rand_seeds = set(random.sample(all_genes, number_of_seed_genes))

        # get rand lcc
        lcc = get_lcc_size(G, rand_seeds)
        l_list.append(lcc)
            

    # get the actual value
    lcc_observed = get_lcc_size(G,gene_set)

    # get the lcc z-score:
    l_mean = np.mean(l_list)
    l_std  = np.std(l_list)
    z_score = (1.*lcc_observed - l_mean)/l_std
    pval = 2*(1-stats.norm.cdf(abs(z_score)))
    

    results = {'LCC_list': l_list,
               'mean': l_mean,
               'std': l_std,
               'z_score': z_score, 
               'p_value':pval}

    return results
