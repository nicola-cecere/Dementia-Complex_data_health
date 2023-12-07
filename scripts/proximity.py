#! /usr/bin/env python

"""
# -----------------------------------------------------------------------
#
# proximity.py
#
# by Emre Guney
# ADAPTED VERSION (2023-11-06)
# 
# -----------------------------------------------------------------------
"""


import networkx as nx
import numpy as np
from scipy import stats

import localization



# =============================================================================
def calculate_closest_distance(net,A,B):
    
    values_outer = []
    
    for node_from in A:
        values = []
        for node_to in B:
            if node_from != node_to:
                try:
                    val = nx.shortest_path_length(net, node_from, node_to)
                    values.append(val)
                except nx.NetworkXNoPath:
                    continue
            else:
                val = 0
                values.append(val)

        if values:
            d = min(values)
            values_outer.append(d)
        
    if len(values_outer) > 0:
        d = np.mean(values_outer)
    else:
        d=''

    return d
    
# =============================================================================
def get_proximity(G, genes, targets, sims):    

    # getting all genes in the network  
    all_genes = sorted(G.nodes())
    gene_set = set(genes) & set(G.nodes())
    target_set = set(targets) & set(G.nodes())
    
    pr_obs = calculate_closest_distance(G, gene_set, targets) 
    if pr_obs == '':
        print('No path between the drug and disease module')
        return None

    
    pr_list  = []

    bucket = localization.group_nodes_by_degree(G)
    # simulations with randomly distributed seed nodes
    for i in range(sims):
        rand_seeds = localization.sample_preserving_degrees(G, gene_set, bucket)
        # get rand proximity
        pr = calculate_closest_distance(G, rand_seeds, target_set)
        pr_list.append(pr)


    # get the z-score:
    pr_mean = np.mean(pr_list)
    pr_std  = np.std(pr_list)
    z_score = (pr_obs - pr_mean)/pr_std
    pval = 2*(1-stats.norm.cdf(abs(z_score)))
    

    results = {'proximity': pr_obs,
               'proximity_list': pr_list,
               'mean': pr_mean,
               'std': pr_std,
               'z_score': z_score, 
               'p_value':pval}

    return results
