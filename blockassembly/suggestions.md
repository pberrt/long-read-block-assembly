# Possibles upgrades

## General

- Seprate the workflow in parts, visually and functionnaly (function of a part in a different file)
- Combine different scripts according to the functions, without loosing previous work.
- Build a SnakeMake file to combine Python scripts and MSA/minimap.
- Define the parameters for benchmark (cliping, multi-k, gene-margin...)
- Define the metrics for benchmark, especially time and RAM consumption. But also k-min to prevent local cycles
- Test on Zamin data

## Input Data

- Value annotation gap by using *None* block.

## Connected components

- At the moment, only the biggest principle component is processed. Need to find a way to apply on each connected component when it is detatched from the block.
- Need high k to detatch the components. Can we cluster the reads based on the kmers (small ones like k=1-3) and then apply the method on each cluster ?

## DBG construction
- DBG construction by iterating on kmers and only read-coherent edges (TwoPaco)
- Compation algorithm in *LaTex*.

## Graph processing

- multi-k improves read coverage (cyclicity of the chromosome can be kept for example) but requires to prune the tips to have longer unitigs
-some kmer of abundance 1 prevent cyclicity, and thus requires to highly increase k if not removed
-remove tips on length and abundance can improve cyclicity and the length of the unitigs.

## Cylcicity
- Self-revcomp-looping tips are not removed and count as a cycle. It occurs with the presence of palindroms (Sample 57) 
- When only 1 *main* cycle pruning is also required to create long unitigs in the main path and thus making it possible to cut at the right place. This can also be improve by checking the cycle itself and identify the *stables* unitigs as it is done later in the workflow
- Improve cycle finding algorithm

## DAG

- Check how to align reads on DAG according to the POA method
- DAG graph are not saved, preventing to use its information for the final consensus

## MSA

- See effect of gene-margin on consensus

## Consensus

- Use the structure of the DAG to reconstruct, not only the topological order
- Write a script for consensus


## Results analysis

- Ability to compute if a reference sequence is included in the *cyclic dag*, or even another *cyclic dag*.
- Compare MR genes fasta on final consensus
- BWA reads on alignement to check where there are mismatch in the dag




