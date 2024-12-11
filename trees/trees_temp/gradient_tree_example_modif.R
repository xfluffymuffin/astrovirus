library('ggtree')
library(ggplot2)
library('phytools')
source("C:\\Users\\gdvov\\OneDrive\\Documents\\gradient_tree_guide\\add_gradient_colors.R")

setwd("C:\\Users\\gdvov\\OneDrive\\Documents\\GitHub\\astrovirus-main\\trees\\iqtree_output\\tempest_analysis\\all_serotypes\\all_serotypes_4948-5720")

# read tree from file
tree = read.tree("all_serotypes_4948-5720.fas.treefile")
tree_rooted = midpoint.root(tree)

# read metadata. info is a dataframe
info = read.csv("metadata.csv")
info = info %>% arrange(GBAC)
write.csv(info, "metadata_sorted.csv", row.names = F)
# run this line to print table




# create basic tree object
t = ggtree(tree_rooted, size=0.75) # size is for lines thickness
t # to see the plot, run this line

# now we will add tip labels
t = ggtree(tree_rooted, size=0.75) + geom_tiplab(size =1)

# show internal node numbers
ggtree(tree_rooted) + geom_text(aes(label=node), hjust=-.3)

# lets color tip labels by host annotation from info dataframe (host column from dataframe)
t = ggtree(tree_rooted, size=0.75) %<+% info + geom_tiplab(size = 2, aes(color=host))
t 

# now let's color the taxa labels by gradient ("code" column in info table defines the order, "color_orf1b16" defines the manual colors)
# but we will have to remove the legend
t = ggtree(tree_rooted, size=0.75) %<+% info + geom_tiplab(size =1, aes(color=label)) +
  scale_color_manual(values=info$color_orf1b16) + geom_treescale() + theme(legend.position = "none")
t


# get taxa names for the whole tree
taxa_names = get_taxa_name(t)
# save taxa names to file
write.table(taxa_names,file="tree_order.csv",col.names=F, row.names=F)

write.table(get_taxa_name(t, 122),file="tree_order_122.csv",col.names=F, row.names=F) # 1
write.table(get_taxa_name(t, 85),file="tree_order_85.csv",col.names=F, row.names=F) # 4 (+ the only seq from 8)
# write.table(get_taxa_name(t, 151),file="tree_order_151.csv",col.names=F, row.names=F) # 8 (+ one seq from 4)
write.table(get_taxa_name(t, 140),file="tree_order_140.csv",col.names=F, row.names=F) # 2
write.table(get_taxa_name(t, 100),file="tree_order_100.csv",col.names=F, row.names=F) # 3
write.table(get_taxa_name(t, 99),file="tree_order_99.csv",col.names=F, row.names=F) # 6
write.table(get_taxa_name(t, 91),file="tree_order_91.csv",col.names=F, row.names=F) # 5


#metadata with gradient colors for order_files.txt (contains "tree_order_657.csv" and "tree_order_811.csv")
info_upd = add_colors2meta("order_files.txt", "metadata_sorted.csv")
write.csv(info_upd, "metadata_upd.csv",row.names=F)

# visualize gradient coloring for two clades
ggtree(tree_rooted, size=0.75) %<+% info_upd + geom_tiplab(size =1, aes(color=label)) +
  scale_color_manual(values=info_upd$color) + geom_treescale() + theme(legend.position = "none")
ggsave("tree_plot.png", plot = last_plot(), width = 110, height = 10, units = "cm", dpi = 600)



# If you need to visualize lots of trees in the same way, you can write a function

#' Function to plot the phylogenetic tree
#' @param tree_file  path to tree file in nwk format
#' @param meta  path to table with metadata
#' @return ggtree object

plot_tree = function(tree_file, meta){
  
  tree =  read.tree(tree_file)
  # root by midpoint
  tree_rooted = midpoint.root(tree)
  # read csv file with metadata
  info = read.csv(meta)
  
  t = ggtree(tree_rooted, size=0.75) %<+% info + geom_tiplab(size =1, aes(color=label)) +
    scale_color_manual(values=info$color) + geom_treescale() + theme(legend.position = "none")
  
  return(t)
  
}

# Then you can get the list of trees in a specific directory (working directory in this example)
trees = list.files(path = "C:\\Users\\gdvov\\OneDrive\\Documents\\gradient_tree_guide\\all_other_trees", full.names = TRUE, pattern = ".treefile")
meta_path = "metadata_upd.csv"
# Create visualization in a loop
for (file in trees){
  # plot tree without legend, add title to plot
  g = plot_tree(file,
                meta_path)  +
    ggtitle(strsplit(strsplit(file, '/')[[1]][2], '.nwk')[[1]][1])
  # save figure in png format. svg, pdf are also possible
  ggsave(paste(file, "png", sep="."),height=15, width=110, units='cm', dpi=1200)
}




