# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 15:32:53 2023


@author: ecsLab

"""

import numpy as np
import pandas as pd


geneBaseN = 4   # motion types
geneLength = 16 

def checkState(i, j, world, worldSize):
    ''' check the circustance of specify position.
        the map is designed such that the left border coincides with the right border.
        and the upper border coincides with the lower border.
    '''
    if i == -1:
        i = worldSize - 1
    elif i == worldSize:
        i = 0
    if j == -1:
        j = worldSize - 1
    elif j == worldSize:
        j = 0
    return world[i][j]


def mutate1spot(gene):
    ''' change 1 gene (randomly seleted) on the chromosome
    '''
    pos = np.random.randint(geneLength)
    gene[pos] = np.random.randint(geneBaseN)
    return gene, pos


def crossover(g1, g2):
    ''' cross-over at a randomly selected point
    '''
    pos = np.random.randint(1, geneLength)
    g1new = np.append(g1[:pos], g2[pos:])
    g2new = np.append(g2[:pos], g1[pos:])
    return g1new, g2new, pos


def geneScore(gene, world, worldSize, Nstep, record=False):
    ''' traversing a chromosome on the world(map) at specified steps(Nstep)
    '''
    posi, posj = 0, 0
    world[posi][posj] = 1
    posiList, posjList = [posi], [posj]
    action_N = []
    for step in range(1, Nstep):
        North  = checkState(posi-1, posj,   world, worldSize)
        East   = checkState(posi,   posj+1, world, worldSize)
        South  = checkState(posi+1, posj,   world, worldSize)
        West   = checkState(posi,   posj-1, world, worldSize)
        
        loc = North*8 + East*4 + South*2 + West
        action = gene[loc]
        
        action_N.append(loc)
        
        if action == 0:   # go north 
            posi -= 1
        if action == 1:   # go east
            posj += 1
        elif action == 2: # go south
            posi += 1
        elif action == 3: # go west
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
        if record:
            posiList.append(posi)
            posjList.append(posj)
    if record:
        return sum(sum(world)), posiList, posjList
    else:
        return sum(sum(world)), action_N

def fitness_scaling(score_all_chro):
    ''' normalizing the fitness of the chromesomes
    '''
    score_all_chro = pd.Series(score_all_chro)
    score_all_chro_rank = score_all_chro.rank(ascending=False, method='dense').values
    score_all_chro_sca = [1/np.sqrt(i) for i in score_all_chro_rank]
    score_nor = score_all_chro_sca/sum(score_all_chro_sca) 
    
    return score_nor


# if __name__ == "__main__":
    
