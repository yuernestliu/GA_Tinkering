# -*- coding: utf-8 -*-
"""

Traversing the whole solution space to find all chromesomes that could have full scores

The output of the script is a txt file with all full-score chromesomes recorded.


@author: ecsLab
"""

import numpy as np
import copy
from math import pow
from base import checkState


gene_num = 4   # the kinds of actions
worldSize = 5  # size of the world(map)

def geneScore(gene, world, worldSize, Nstep): 
    posi, posj = 0, 0
    world[posi][posj] = 1
    action_N = []
    for step in range(1, Nstep):
        North  = checkState(posi-1, posj,   world, worldSize)
        East   = checkState(posi,   posj+1, world, worldSize)
        South  = checkState(posi+1, posj,   world, worldSize)
        West   = checkState(posi,   posj-1, world, worldSize)
        
        action_N.append(North*8 + East*4 + South*2 + West)
        action = gene[North*8 + East*4 + South*2 + West]
        
        if action == 0:
            posi -= 1
        if action == 1:
            posj += 1
        elif action == 2:
            posi += 1
        elif action == 3:
            posj -= 1

        if posi == -1:
            posi = worldSize - 1
        elif posi == worldSize:
            posi = 0
        if posj == -1:
            posj = worldSize - 1
        elif posj == worldSize:
            posj = 0
        world[posi][posj] = 1
    return sum(sum(world)), action_N


def evaluate_each_chrom(filename):
    '''transverse the whole solution space and record the chromosomes with full-score
    filename: txt file name
    '''
    NstepRobot = worldSize*worldSize  # total steps of the roboat
    world0 = np.zeros((worldSize, worldSize), dtype='int8')
    
    chrom_num = 0  # count the number of chromesomes that have full score
    # chrom_len_all = []  # record the effective length of each full-score chromosome
    for l0 in range(gene_num):
        
        for l1 in range(gene_num):
            
            for l2 in range(gene_num):
                
                for l3 in range(gene_num):
                        
                    for l4 in range(gene_num):
                        
                        for l5 in range(gene_num):
                            
                            for l6 in range(gene_num):
                                
                                for l7 in range(gene_num):
                                    
                                    for l8 in range(gene_num):
                                        
                                        for l9 in range(gene_num):
                                            
                                            for l10 in range(gene_num):
                                                
                                                for l11 in range(gene_num):
                                                    
                                                    for l12 in range(gene_num):
                                                        
                                                        for l13 in range(gene_num):
                                                            
                                                            for l14 in range(gene_num):
                                                                
                                                                for l15 in range(gene_num):
                                    
                                                                    chrom = [l0,l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,l14,l15]
                                                                    score, _ = geneScore(chrom, copy.deepcopy(world0), worldSize, NstepRobot)
                                                                        
                                                                    if score == NstepRobot:
                                                                        chrom_num+=1
                                                                        # chrom_len_all.append(len(np.unique(action_N).tolist()))
            
                                                                        with open(file=filename, mode='a',encoding='utf-8') as f:
                                                                            f.write(str(chrom)+'\n')
                                                                    
    return chrom_num


if __name__ == "__main__":
    filename = r'D:\chroms.txt'
    chrom_num = evaluate_each_chrom(filename)
    print(f'full score chromesomes: {chrom_num} ')
    # save chrom_len_all