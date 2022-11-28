# InverseMappingFunctions

This repository contains the data and scripts to reproduce the inverse mapping function analyses presented in Kivikoski *et al.* (2022) *Predicting recombination frequency from map distance*.

## Structure of the repository

### Analyses

Analyses steps as R scripts.

### Data

- Nine-spined stickleback (*Pungitius pungitius*) Marey map and offspring haplotypes.
- Three-spined stickleback (*Gasterosteus aculeatus*) Marey map and offspring haplotypes.

### Examples

 Small example data and a script to demonstrate how to infer the likelihoods for different bivalent crossover counts and how to predict the recombination frequency between two markers, by using the *p<sub>0</sub>(k)* function or inverse mapping functions.

 To run the example, type at the command prompt: 

 `git clone https://github.com/mikkokivikoski/InverseMappingFunctions.git
  
  cd InverseMappingFunctions/Examples/
  
  Rscript example.r`

### Functions

- emAlgortihmForCoProbabilities.r

  R implementation of the expectation-maximization algorithm to infer bivalent crossover frequencies from gametic recombinations (see Yu, K. and Feingold, E. (2001), Estimating the frequency distribution of crossovers during meiosis from recombination data. *Biometrics*, 57: 427-434.).

- p0kFunction.r 

  R implementation of the *p<sub>0</sub>(k)* function developed in the manuscript.
