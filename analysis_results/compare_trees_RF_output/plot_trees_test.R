library(ggtree)
library(treeio) #read.beast
library(ape) # for phylo objects
library(tidyr)

# set  working directory
setwd("C:\\Users\\gdvov\\OneDrive\\Documents\\GitHub\\astrovirus-main\\analysis_results\\compare_trees_RF_output")
# Load functions from source code
source("get_subtrees.R")

# converts file with common subtrees in newick format into list with taxa names
# *_commontrees.txt" is a file with common subtrees between two phylogenetic trees
# Each subtree is in newick format
# this file is output of get_RF_halflife.py
# Example: python.exe .\get_RF_halflife.py -tree1 .\norovirus_vp1.tree -tree2 .\norovirus_rdrp_g2.tree

subtrees = get_subtrees("serotype_1\\5869-7192_vs_ORF2_ml\\serotype_1_5869-7192_strict_skyline_400_combined_ORF_2_cleared_alignment_RNA_serotypes_no_amb.fas_commontrees.txt")

# Path to file with tree in nexus format produced by BEAST
tree1_name = "input_trees\\Bayesian_regional\\serotype_1_5869-7192_strict_skyline_400_combined.tree"
tree1 = read.beast(tree1_name)

# Plot beast subtree and color braches of common subtrees
p <- ggtree(tree1,mrsd="2025-01-01")
groupOTU(p, subtrees, 'Clade') + aes(color=Clade) +
  theme(legend.position="right") + scale_color_manual(values=c("black", "firebrick"))

