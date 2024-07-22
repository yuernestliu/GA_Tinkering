# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 17:23:08 2023

gengerate the evolution tree of the full-score chromosome

The 'evolution' function below will store the evolution information,
and then it would be visualized by 'GA_graph' function

@author: ecslab

"""
import numpy as np
import graphviz
import os
import time
import copy
import random
from base import checkState, mutate1spot, crossover, geneScore, fitness_scaling


geneBaseN = 4 
geneLength = 16 

def GA_graph(target_chr_idx, geneComp_all, scoreArray_all):
    '''
    Input：
        1. target_chr_idx：the index of full-score chromosome (in the chromosome pool)
        2. geneComp_all: multi-level list included the evolution information  [[father1_id, father2_id，cross_loc],[muta_father_id，muta_loc], father_id] 
        3. scoreArray_all： array size = generation_times * gene_pool_size
    
    About the graph：
    node：the name of the node use the generation index and the chromosome index. e.g. 1_2 represents the 2nd chromosome in the 1st generation
    edge: represents the relationship between the chromosome, and 
        1. black：the same
        2. green：mutation
        3. red：cross over
        4. blue：cross over & mutation
    
    color of the node：gray scale, the deeper, the higher score
    
    return: dot
    '''

    dot = graphviz.Digraph('evolution-tree') 
    # generation index start from 0
    totol_gen = len(geneComp_all)                        # totoal evolution times
    score_max = max([max(i) for i in scoreArray_all])    # maximum score
    final_node = str(totol_gen)+'_'+str(target_chr_idx)
    score = scoreArray_all[totol_gen, target_chr_idx]
    color = '0 0 '+str(1-score/(2*score_max))
    dot.node(final_node, style='filled', color=color)    # node name
    
    edge_type = [0,0,0,0]  # same, mutation, cross-over, cross&mutation
    for gen in np.arange(totol_gen, 0, -1, dtype=int):   # generation index: start from totol_gen, step by -1
        
        if gen == totol_gen:  # start from target chromosome
            c_idx = [target_chr_idx]                     # chromosome index in the current generation
            last_nodes = [final_node]                    # node name of the former generation
        else:
            # delete the repeat chromosomes
            new_idx_save_uniq = []
            last_nodes_save_uniq = []
            for i in range(len(new_idx_save)):
                if new_idx_save[i] not in new_idx_save_uniq:
                    new_idx_save_uniq.append(new_idx_save[i])
                    last_nodes_save_uniq.append(last_nodes_save[i])    
                
            c_idx = new_idx_save_uniq
            last_nodes = last_nodes_save_uniq
        
        new_idx_save = []  # add new node
        last_nodes_save = []
        
        for idx, last_node in zip(c_idx, last_nodes):        # when multiple node exist in one generation
        
            gen_info = geneComp_all[gen-1][idx]   # 
                
            if gen_info[0] is None and gen_info[1] is None:
                edge_type[0]+=1
                last_idx = gen_info[2]
                new_idx_save.append(last_idx)
                new_node = str(gen-1)+'_'+str(last_idx)
                last_nodes_save.append(new_node)
                
                score = scoreArray_all[gen-1, last_idx]
                color = '0 0 '+str(1-score/(2*score_max))    # "0.000 0.000 0.000"	"black"  "0.000 0.000 1.000"	"white"
                
                dot.node(new_node, style='filled', color=color)
                dot.edge(new_node, last_node, color='black')
                
            elif gen_info[0] is None and type(gen_info[1]) is list:  # mutation only
                edge_type[1]+=1
                # last_idx = chosenMap[gen-1][idx]
                last_idx = gen_info[1][0]
                new_idx_save.append(last_idx)
                new_node = str(gen-1)+'_'+str(last_idx)
                last_nodes_save.append(new_node)
                
                score = scoreArray_all[gen-1, last_idx]
                color = '0 0 '+str(1-score/(2*score_max))
                
                dot.node(new_node, style='filled', color=color)
                dot.edge(new_node, last_node, color='green')
            
            elif type(gen_info[0]) is list and gen_info[1] is None:  # cross-over only
                edge_type[2]+=1
                # father 1
                last_idx1_1 = gen_info[0][0]
                new_idx_save.append(last_idx1_1)
                new_node1 = str(gen-1)+'_'+str(last_idx1_1)
                last_nodes_save.append(new_node1)
                
                score = scoreArray_all[gen-1, last_idx1_1]
                color = '0 0 '+str(1-score/(2*score_max))
                
                dot.node(new_node1, style='filled', color=color)
                dot.edge(new_node1, last_node, color='red')
                
                # father 2
                last_idx2_1 = gen_info[0][1]
                new_idx_save.append(last_idx2_1)
                new_node2 = str(gen-1)+'_'+str(last_idx2_1)
                last_nodes_save.append(new_node2)
                
                score = scoreArray_all[gen-1, last_idx2_1]
                color = '0 0 '+str(1-score/(2*score_max))
                
                dot.node(new_node2, style='filled', color=color)
                dot.edge(new_node2, last_node, color='red')
            
            else:  # cross-over and mutaion
                edge_type[3]+=1
                # father 1
                last_idx1_1 = gen_info[0][0]
                new_idx_save.append(last_idx1_1)
                new_node1 = str(gen-1)+'_'+str(last_idx1_1)
                last_nodes_save.append(new_node1)
                
                score = scoreArray_all[gen-1,last_idx1_1]
                color = '0 0 '+str(1-score/(2*score_max))
                
                dot.node(new_node1, style='filled', color=color)
                if gen_info[1][0] == gen_info[0][0]:
                    dot.edge(new_node1, last_node, color='blue')
                else:
                    dot.edge(new_node1, last_node, color='red')
                    
                # father 2
                last_idx2_1 = gen_info[0][1]
                new_idx_save.append(last_idx2_1)
                new_node2 = str(gen-1)+'_'+str(last_idx2_1)
                last_nodes_save.append(new_node2)
                
                score = scoreArray_all[gen-1, last_idx2_1]
                color = '0 0 '+str(1-score/(2*score_max))
                
                dot.node(new_node2, style='filled', color=color)
                if gen_info[1][0] == gen_info[0][1]:
                    dot.edge(new_node2, last_node, color='blue')
                else:
                    dot.edge(new_node2, last_node, color='red')
    return dot,edge_type


def evolution(crossoverProb, mutationProb):
    '''evolution with specified cross-over and mutation rate
    
    Output：
        1. geneComp_all: evolution information, [[index of father1, index of father2，cross-over position],[index of father chromosome，mutation position], index of father chromosome] 
        2. scoreArray_all： array - generation_times * gene_pool_size
        3. best_idx: index of the full-score chromosomes when evolution stop
        4. best_length: the effective length of the full-score chromosome
    '''
    
    Ngeneration = 200  
    genePoolSize = 20    
    
    worldSize = 5 
    NstepRobot = worldSize*worldSize  

    genePoolIndex = [i for i in range(genePoolSize)]
    world0 = np.zeros((worldSize, worldSize), dtype='int8')
    genePool = []
    # scoreArray_all = np.zeros((Ngeneration, genePoolSize), dtype=float)
    scoreArray_all = []
    
    # 第0代
    for i in range(genePoolSize):
        genePool.append(np.random.randint(low=0, high=geneBaseN, size=geneLength, dtype='int8'))
    scoreArray = np.array( [0.0]*genePoolSize )
    locArray = [[] for _ in range(genePoolSize)]
    
    for i, gene in enumerate(genePool):
        scoreArray[i], locArray[i] = geneScore(gene, copy.deepcopy(world0), worldSize, NstepRobot)  # copy.deepcopy(world0)
    scoreArray_all.append(copy.deepcopy(scoreArray))
    
    # gengeration 1 to end
    geneComp_all = []
    
    for generation in range(1, Ngeneration+1):
        gain = 1.5
        use_gene_num = [len(np.unique(locArray[i]))*gain for i in range(len(scoreArray))]
        scoreArray_new = scoreArray + use_gene_num
        
        scoreArray_norm = fitness_scaling(scoreArray_new)
        indexChosen = random.choices(genePoolIndex, weights = scoreArray_norm, k = genePoolSize)
        genePoolNew = []
        for i in indexChosen:
            genePoolNew.append(genePool[i].copy())
    
        # Cross over and mutation
        geneComp =  [[None,None,indexChosen[i]] for i in range(genePoolSize)] 
        for i in range(0, genePoolSize, 2):
            if random.random() < crossoverProb:
                genePoolNew[i], genePoolNew[i+1], pos = crossover(genePoolNew[i], genePoolNew[i+1])
                geneComp[i][0] = [indexChosen[i], indexChosen[i+1], pos]
                geneComp[i+1][0] = [indexChosen[i+1], indexChosen[i], pos]
                
            if random.random() < mutationProb:
                genePoolNew[i], pos = mutate1spot(genePoolNew[i])
                geneComp[i][1] = [indexChosen[i], pos]
                
            if random.random() < mutationProb:
                genePoolNew[i+1], pos = mutate1spot(genePoolNew[i+1])
                geneComp[i+1][1] = [indexChosen[i+1], pos]
                
        genePool = genePoolNew
        # geneComp_all[generation-1] = geneComp
        geneComp_all.append(geneComp)
        
        for i, gene in enumerate(genePool):
            scoreArray[i], locArray[i] = geneScore(gene, copy.deepcopy(world0), worldSize, NstepRobot)
        # scoreArray_all[generation,:] = scoreArray
        scoreArray_all.append(copy.deepcopy(scoreArray))
        
        # stop when full-score chromosome exist
        best_idx = []
        if np.max(scoreArray) == worldSize*worldSize:
            evo_times = generation
            best_idx = [i for i in range(len(scoreArray)) if scoreArray[i]==max(scoreArray)]  # store all full-score chromosome
            break
    
    all_len = [len(np.unique(i)) for i in locArray]
    
    return geneComp_all, np.array(scoreArray_all), best_idx, [all_len[i] for i in best_idx]


if __name__ == "__main__":
    crossoverProb = 0.2
    mutationProb = 0.2
    
    run_times = 1
    times = 0
    while times < run_times:
        times +=1
        geneComp_all, scoreArray_all, best_idx, best_len = evolution(crossoverProb, mutationProb)
        
        if len(best_idx)>0:
            print(best_idx)
            dot, edge_type = GA_graph(best_idx[0], geneComp_all, scoreArray_all)  # notice that plot the first one as default
            
            #%%  
            dot.format = 'png'
            save_path = r'Trees\\'+str(best_len[0])  # make new folder under current dir
            if not os.path.exists(save_path):
                os.makedirs(save_path)
            
            dot.render(directory=save_path).replace('\\', '/')
