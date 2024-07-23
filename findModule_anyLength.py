# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 12:49:08 2024

Find modules with different number of genes (Figure 6.)
Different number of genes in m could be set at line 57 (e.g. 2,3,4)

@author: ecslab
"""

import copy
import numpy as np
import random
import matplotlib.pyplot as plt
# import sys
# sys.path.append(r'D:\Nutstore\Project code\GA_v2\Span_1')  # 
from Robot_span_1 import geneScore

#%% Basic parameters
worldSize = 10
NstepRobot = worldSize*worldSize
world0 = np.zeros((worldSize, worldSize), dtype='int8')
geneBaseN = 4   # types of motions
geneLength = 16 
simu_times = 1000  # equals to the size of population
np.random.seed(0)
chrome_pool = np.random.randint(low=0, high=geneBaseN, size=(simu_times,geneLength), dtype='int8')

#%% scores of the original chromosomes
score_ori = []
for chrome in chrome_pool:
    world = copy.deepcopy(world0)
    score, _ = geneScore(chrome, world, worldSize, NstepRobot)
    score_ori.append(score)

#%%    
movement = {
    0: [0,1,2,3],
    1: [0,1,2],
    2: [0,1,3],
    3: [0,1],
    4: [0,2,3],
    5: [0,2],
    6: [0,3],
    7: [0],
    8: [1,2,3],
    9: [1,2],
    10: [1,3],
    11: [1],
    12: [2,3],
    13: [2],
    14: [3]
    #     15: []
}

pre_module_pool_size = 1000  # number of the random module
gene_num = 2 # how many genes in a module

# generate gene locations
gene_locs = []  # dimension: pre_module_pool_size * gene_num
for i in range(pre_module_pool_size):
    np.random.seed(i)
    gene_locs.append(np.random.choice(geneLength-1, gene_num, replace=False).tolist())
#     if len(gene_locs) > 1:
#         for j in gene_locs:
#             while set(np.random.choice(geneLength-1, gene_num, replace=False).tolist()) == set(j):
#                 np.random.seed(i+pre_module_pool_size)
#                 break
#         gene_locs.append(np.random.choice(geneLength-1, gene_num, replace=False).tolist())
#     else:
#         gene_locs.append(np.random.choice(geneLength-1, gene_num, replace=False).tolist())
gene_locs = np.array(gene_locs)

# Basic units of a module: [[location1, motion1],[location2, motion2]...]
# Generate the pool of modules
module_ini = [[0,0] for _ in range(gene_num)]  # [gene_location, movement]
pre_module_pool = [module_ini for _ in range(pre_module_pool_size)]
pre_module_pool = np.array(pre_module_pool)

set_seed = 0
for i,pre_module in enumerate(pre_module_pool):
    
    for ii,gene_loc in enumerate(gene_locs[i]):
        set_seed += 1
        pre_module[ii,0] = gene_loc
        movement_range = movement[gene_loc]
        random.seed(set_seed)
        pre_module[ii,1] = random.choice(movement_range)

# compute the score after module replacement (1. replace all genes at one time 2. replace 1 gene of a module at a time)
pre_module_scores = np.zeros((pre_module_pool_size, simu_times), dtype=int)
sub_scores = np.zeros((pre_module_pool_size, simu_times, gene_num), dtype=int)

for j,pre_module in enumerate(pre_module_pool):
    add_module_score = np.zeros((simu_times,), dtype=int)
    sub_score = np.zeros((simu_times,gene_num), dtype=int)
    
    for i,chrome in enumerate(chrome_pool):
        # replace all
        chrome_1 = copy.deepcopy(chrome)
        for ii in range(gene_num):
            chrome_1[pre_module[ii][0]] = pre_module[ii][1]
        score, _ = geneScore(chrome_1, copy.deepcopy(world0), worldSize, NstepRobot)
        add_module_score[i] = score
        
        # replace 1 gene at a time
        chrome_2 = copy.deepcopy(chrome)
        for jj in range(gene_num):
            chrome_2[pre_module[jj][0]] = pre_module[jj][1]
            score, _ = geneScore(chrome_2, copy.deepcopy(world0), worldSize, NstepRobot)
            sub_score[i,jj] = score

    pre_module_scores[j,:] = add_module_score
    sub_scores[j,:,:] = sub_score

# Calculate the improvement in total score when replacing multiple genes at the same time, 
# compared to the replacement of individual genes, as an evaluation for the module   
score_ori_mean = np.mean(score_ori)  # 
gain_module = np.mean(pre_module_scores,1) - score_ori_mean
sub_gain_sum = 0
for i in range(gene_num):
    sub_gain_sum += np.mean(sub_scores[:,:,i],1) - score_ori_mean
    
pre_module_score = gain_module - sub_gain_sum


#%% Plot the result
fig = plt.figure()

# Compute histogram without plotting
counts, bins = np.histogram(pre_module_score, bins=30)

# Calculate relative frequencies
relative_frequencies = counts / counts.sum()

# Plot histogram using relative frequencies
plt.bar(bins[:-1], relative_frequencies, width=np.diff(bins), edgecolor='black', align='edge', color='royalblue')
plt.title("The length of $m$ = 2", fontsize=17)
plt.xlabel("Modular index $\mu$", fontsize=15)
plt.ylabel('Proportion', fontsize=15)

plt.savefig(r'm=2.png', dpi=300)

# Show the plot
plt.show()