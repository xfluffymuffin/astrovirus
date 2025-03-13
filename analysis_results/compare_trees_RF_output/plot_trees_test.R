library(ggtree)
library(treeio) #read.beast
library(ape) # for phylo objects
library(tidyr)
# Load functions from source code
source("get_subtrees.R")

# converts file with common subtrees in newick format into list with taxa names
# *_commontrees.txt" is a file with common subtrees between two phylogenetic trees
# Each subtree is in newick format
# this file is output of get_RF_halflife.py
# Example: python.exe .\get_RF_halflife.py -tree1 .\norovirus_vp1.tree -tree2 .\norovirus_rdrp_g2.tree

subtrees = get_subtrees("1330_2896_2908_4090_commontrees.txt")

# Path to file with tree in nexus format produced by BEAST
tree1_name = "1330_2896.tree"
tree1 = read.beast(tree1_name)

# Plot beast subtree and color braches of common subtrees
p <- ggtree(tree1,mrsd="2025-01-01")
groupOTU(p, subtrees, 'Clade') + aes(color=Clade) +
  theme(legend.position="right") + scale_color_manual(values=c("black", "firebrick"))

