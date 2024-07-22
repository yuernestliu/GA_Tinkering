# Out of Randomness: How Evolution Benefits from Modularity

## Introduction

This repository contains the source code accompanying our paper "Out of Randomness: How Evolution Benefits from Modularity". In this study, we explore how evolution-like processes, particularly those involving mutation and crossover, can effectively search for complex solutions.

## Project Purpose

The main objectives of this project are:

1. To quantitatively compare the effectiveness of Genetic Algorithms (GAs) against brute force random searches in finding solutions, especially complex ones.

2. To investigate potential biases and asymmetries in GA-based searches, particularly towards more intricate solutions.

3. To examine how the grouping of module components affects the formation of larger, more complex solutions.

4. To draw parallels between our findings and patterns observed in real biological systems, such as the clustering of gene-encoding sequences in prokaryotic organisms.

## Key Findings

Our experiments reveal that:

- GAs significantly improve the probability of finding solutions compared to random searches.
- There's a bias towards more complex solutions in GAs, likely due to the crossover process.
- Grouping module components aids in forming larger, more complex solutions.

These findings highlight the importance of modularity and modular recombination in enhancing solution space exploration, both in evolutionary algorithms and potentially in natural evolution.

## Repository Contents
- `base.py`: Contains basic functions for chromosome manipulation, world state checking, and fitness calculation.
- `all_full_score_solutions.py`: Script to find all possible full-score solutions by exhaustive search.
- `evolution_tree.py`: Main script for running the genetic algorithm and generating evolution trees (for Fig. S3 in Supplementary Information of the paper).
- `evolution_tree_prune.py`: Pruning the extra edges from the graphs (for Fig.5a in the paper).


## Requirements
This project requires Python 3.x and the following Python libraries:
- NumPy
- graphviz
