# compare_trees_RF
compare phylogenetic trees and calculate TMRCAs of recombinant forms 

Functions for identifing and visualizing common subtrees in two phylogenetic trees. One of the trees should be a time-scaled phylogeny produced by BEAST.


## get_RF_halflife.py

This script finds common subtrees between two trees and calculates median of their MRCA ages. The subtrees' tMRCA are extracted from the first tree which should be a time-scaled phylogeny produced by BEAST software. The second tree's format could be nexus or newick. Output:
- *_commontrees.txt - text file with common subtrees in newick format
- *_heights.txt - ages of common subtrees

Importantly, the trees' taxa should not necessarily coincide. 

Usage
```
get_RF_halflife.py [-h] -tree1 TREE_BEAST -tree2 TREE2 [-pthr POSTERIOR_THRESHOLD]

options:
  -h, --help            show this help message and exit
  -tree1 TREE_BEAST, --tree_beast TREE_BEAST
                        Path to time tree in nexus format inferred using BEAST software
  -tree2 TREE2, --tree2 TREE2
                        Path to phylogenetic tree in nwk or nexus format
  -method METHOD, --method METHOD
                        Method for determining recombinant forms. 'subtrees', 'bipartitions' or 'all'
  -pthr POSTERIOR_THRESHOLD, --posterior_threshold POSTERIOR_THRESHOLD
                        Threshold for posterior values of nodes to count. Ranges from 0 to 1.
  -bthr BOOTSTRAP_THRESHOLD, --bootstrap_threshold BOOTSTRAP_THRESHOLD
                        Threshold for bootstrap values of branches to count (for 'bipartitions' method). Ranges from 0
                        to 100.
```


## get_subtrees.R

Helper functions for parsing *_commontrees.txt


## plot_trees_test.R

Example of visualization.