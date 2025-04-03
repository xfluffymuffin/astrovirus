import argparse
import hashlib
import re
import sys
import os

import dendropy as dpy
import numpy as np
import pandas as pd

# Regular expression for decimals
locval_re = re.compile(r'[0-9]+\.[0-9E\-]+')

'''
This function traverses beast tree and adds the attributes 'posterior' and 'height' for all tree nodes according to the annotations.

Return tree with updated nodes.
'''

def parse_beast_tree_node_info(tree):

    for nd in tree.postorder_node_iter():
        node_comment = nd._get_annotations()
        for annot in node_comment:
            if str(annot).startswith('posterior='):
                posterior =  float(locval_re.search(str(annot)).group())
                nd.posterior = posterior
            if str(annot).startswith('height='):
                height =  float(locval_re.search(str(annot)).group())
                nd.height = height

    return tree

'''
This function sorts alphabetically taxa names in bipartition produced by split_as_newick_string() (dendropy)

bipart_str - bipartition split represented as a newick string

Returns bipart_str with sorted taxa

Example:
bipart_str = "((B,A),(C));"

sort_bipart(bipart_str)
"((A,B),(C));"
'''
def sort_bipart(bipart_str):
    if '((' in bipart_str:
        bipart_l = bipart_str.split('), (')
        #print(bipart_l)
        for i in range(len(bipart_l)):
            bipart_l[i] = bipart_l[i].strip('((').strip('));')
            if len(bipart_l[i]) == 1:
                continue
            else:
                part_l = bipart_l[i].split(', ')
                part_l.sort()
                bipart_l[i] = ','.join(part_l)
        
        bipart_l.sort()
        bipart_str = '((' + bipart_l[0] + ')' + ',' + '(' + bipart_l[1] + '))'
        return bipart_str
    else:
        return ""

'''
Postorder traverses the tree and adds the attribute 'hash' for each node. 
Hash (sha256) is calculated for string that includes all bipartitions in subtree with the descendants of a particular node.
Returns dictionary with hashes and corresponding subtrees in newick format.
'''

def add_hashes(tree, timescaled=False):
    tree_hashes = {}
    biparts = []
    if timescaled:
        subtree_heights = {}
    for nd in tree.postorder_node_iter():
        if nd.is_leaf():
            continue
        else:
            subtree = dpy.Tree(seed_node=nd.extract_subtree())
            subtree.encode_bipartitions()

            bipart_list = [sort_bipart(bip.split_as_newick_string(subtree.taxon_namespace)) for bip in subtree.encode_bipartitions()]
            bipart_list.sort()
            stri_bipart = ';'.join(bipart_list[1:])

            h = hashlib.new('sha256',usedforsecurity=False)
            #h.update(json.dumps(nd._as_newick_string(edge_lengths=None)).encode('utf-8'))
            h.update((stri_bipart).encode('utf-8'))
            nd.hash = h.hexdigest()
            tree_hashes[nd.hash] = nd._as_newick_string(edge_lengths=None)
            tree_hashes[nd.hash] = re.sub("\)[0-9]+", ")", tree_hashes[nd.hash])
            #if "MK073889_USA_human_2016_GII.4_GII.P4" in nd._as_newick_string(edge_lengths=None):
            #    #print(nd._as_newick_string(edge_lengths=None))
            #    biparts.append(stri_bipart)
            
            if timescaled:
                min_leaf_height = 100
                for leaf in nd.leaf_iter():
                    if leaf.height < min_leaf_height:
                        min_leaf_height = leaf.height
                subtree_heights[nd.hash] = min_leaf_height

    if timescaled:
        return tree_hashes, subtree_heights
    return tree_hashes

'''
Extracts heights of nodes which subtrees correspond in two trees

Input:
tree1 - time inferred using BEAST. The nodes of this tree are supposed to have attribute height and hash. 
The attribute posterior is also needed if the threshold for posterior value is used.

hashes_tree2 - dictionary with hashes for the second tree you would compare with

posterior_thr - threshold posterior value. If posterior of the node is lower that the threshold, the height of subtree will not be counted.

Output:

heights - list with heights of subtrees common between tree1 and tree2

'''
def get_common_subtrees(tree1, subtree_times1, hashes_tree2, posterior_thr = None):
    heights = []
    common_subtrees = []
    stack = [tree1.seed_node]
    while stack:
        node = stack.pop()
        if node.is_leaf():
            continue
        if node.hash in hashes_tree2.keys():
            # check node posterior
            if posterior_thr:
                if node.posterior > posterior_thr:
                    #yield node.height
                    #print(node.height)
                    heights.append(node.height - subtree_times1[node.hash])
                    common_subtrees.append(hashes_tree2[node.hash])
                else:
                    stack.extend(n for n in reversed(node._child_nodes))
            else:
                #print(str(node.height) + '-' + str(subtree_times1[node.hash]))
                heights.append(node.height - subtree_times1[node.hash])
                common_subtrees.append(hashes_tree2[node.hash])
        else:
            stack.extend(n for n in reversed(node._child_nodes))
    return heights, common_subtrees


'''
Splits tree to bipartitions. For newick files, saves bootstrap support of a corresponding bipartitions. For nexus files obtained from BEAST saves posterior support and node height.

Input:
tree - time in newick format or in nexus format (inferred using BEAST).
treetype - type of tree ('newick' or 'nexus')

Output:

bipartitions - dictionary:
    keys - hash of sorted bipartition
    values - list. For newick tree, bipartitions[key]=[bootstrap, sorted bipartition as newick string]
    For nexus tree bipartitions[key]=[posterior, node height, sorted bipartition as newick string]

'''

def encode_bipartitions(tree, treetype="newick"):
    #counter for non binary nodes
    k=0
    bipartitions = {}
    tree.encode_bipartitions(suppress_storage=True)
    for bipart in tree.bipartition_edge_map:
        edge = tree.bipartition_edge_map[bipart]
        #print(edge.head_node)
        #print(bipart.split_as_newick_string(tree.taxon_namespace))
        if edge.is_leaf():
            continue
        else:
            bipart_sorted = sort_bipart(bipart.split_as_newick_string(tree.taxon_namespace))
            h = hashlib.new('sha256',usedforsecurity=False)
            h.update((bipart_sorted).encode('utf-8'))
            if treetype == "newick":
                if edge.head_node.label!= None:
                    support = float(edge.head_node.label)
                else:
                    # Case for non-binary tree, when several leafs cpme from one node
                    support = 0
                    k+=1
    
                #print(bipart.split_as_newick_string(tree.taxon_namespace))
                bipartitions[h.hexdigest()] = [support,bipart_sorted]
            else:
                support = edge.head_node.posterior
                height = edge.head_node.height

                min_leaf_height = 100
                for leaf in edge.head_node.leaf_iter():
                    if leaf.height < min_leaf_height:
                        min_leaf_height = leaf.height
                newick_substr = re.sub("\)[0-9]+", ")", edge.head_node._as_newick_string(edge_lengths=None))  
                bipartitions[h.hexdigest()] = [support, height, height - min_leaf_height,bipart_sorted,newick_substr]
    print("Bipartitions from non-binary nodes = {}".format(k))
    return(bipartitions)


'''
Extracts heights of nodes which bipartitions correspond in two trees

Input:

biparts_tree1 - dictionary for bipartitions of tree1 which is time tree in nexus format inferred using BEAST. 

biparts_tree2 - dictionary for bipartitions of tree2 which can be both in newick and time tree in nexus format

Taxa labels must coincide in two trees!
The times and common subtrees are calculated using tree1.

Output:

heights - list with heights of bipartitions common in tree1 and tree2
subtrees - list with subtrees of common bipartitions in tree1 and tree2

'''


def get_common_biparts(biparts_tree1, biparts_tree2, posterior_thr=None, bootstrap_sup=None):
    heights = []
    subtrees = []
    for h in biparts_tree1:
        if posterior_thr!=None and biparts_tree1[h][0] < posterior_thr:
            continue
        if h in biparts_tree2.keys():
            #print(biparts_tree1[h][0])
            if bootstrap_sup!=None and biparts_tree2[h][0] < bootstrap_sup:
                continue
            #print(biparts_tree1[h][1])
            #print(biparts_tree1[h][2])
            #print(biparts_tree2[h][0])
            heights.append(biparts_tree1[h][2])
            subtrees.append(biparts_tree1[h][4])

    return heights, subtrees





if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-tree1", "--tree_beast", type=str,
                        help="Path to time tree in nexus format inferred using BEAST software", required=True)
    parser.add_argument("-tree2", "--tree2", type=str,
                        help="Path to phylogenetic tree in nwk or nexus format", required=True)
    parser.add_argument("-method", "--method", type=str,
                        help="Method for determining recombinant forms. 'subtrees', 'bipartitions' or 'all' ", required=True)
    parser.add_argument("-pthr", "--posterior_threshold", type=float,
                        help="Threshold for posterior values of nodes to count. Ranges from 0 to 1.")
    parser.add_argument("-bthr", "--bootstrap_threshold", type=float,
                        help="Threshold for bootstrap values of branches to count (for 'bipartitions' method). Ranges from 0 to 100.")

    args = parser.parse_args()

    tree1 = dpy.Tree.get_from_path(args.tree_beast, 'nexus')
    print("Getting annotation of tree nodes...")
    tree1 = parse_beast_tree_node_info(tree1)
    print("Done")
    
    try:
        tree2 = dpy.Tree.get_from_path(args.tree2, 'nexus')
    except:
        try:
            tree2 = dpy.Tree.get_from_path(args.tree2, 'newick')
        except:
            print("Couldn't read tree2!")
            sys.exit(1)
            
    tree1_name = os.path.splitext(os.path.split(args.tree_beast)[-1])[0]
    tree2_name = os.path.splitext(os.path.split(args.tree2)[-1])[0]
    #print(tree1_name)
    
    if args.method == 'bipartitions' or args.method == 'all':
        print("Encoding bipartitions for tree 1...")
        biparts_tree1 = encode_bipartitions(tree1,"nexus")
        print("Done")
        
        print("Encoding bipartitions for tree 2...")
        biparts_tree2 = encode_bipartitions(tree2,"newick")
        print("Done")
        
        heights_bip, subtrees_bip = get_common_biparts(biparts_tree1, biparts_tree2)
        
        df = pd.DataFrame(biparts_tree1).T
        df.columns = ['posterior', 'height_raw', 'height_corrected', 'bip', 'subtree']

        df_nwk = pd.DataFrame(biparts_tree2).T
        df_nwk.columns = ['support', 'bip']
        
        print("Number of bipartitions with posterior > {} = {}".format(str(args.posterior_threshold),df[df['posterior']>args.posterior_threshold].shape[0]))
        print("Number of bipartitions with bootstrap > {} = {}".format(str(args.bootstrap_threshold), df_nwk[df_nwk['support']>args.bootstrap_threshold].shape[0]))
        
        print("Number of coinciding bipartitions in tree1 and tree2 with no thresholds {}".format(len(subtrees_bip)))
        print("The median height of common bipartitions with no thresholds is {}".format(round(np.median(heights_bip),4)))


        heights_bip, subtrees_bip = get_common_biparts(biparts_tree1, biparts_tree2, args.posterior_threshold, args.bootstrap_threshold)        
        print("Number of coinciding bipartitions in tree1 and tree2 {}".format(len(subtrees_bip)))
        print("The median height of common bipartitions is {}".format(round(np.median(heights_bip),4)))
    
    
    
        with open(os.path.join(os.getcwd(),tree1_name + '_' + tree2_name + "_commontrees_bip.txt"), 'w') as file:
            file.write("\n".join(subtrees_bip) + "\n")
        file.close()

        with open(os.path.join(os.getcwd(), tree1_name + '_' + tree2_name + "_heights_bip.txt"), 'w') as file:
            file.write("\n".join([str(h) for h in heights_bip]) + "\n")
        file.close()

    if args.method == 'subtrees' or args.method == 'all':
        print("Calculating hashes for tree 1...")
        hashes_tree1, subtree_times1 = add_hashes(tree1, timescaled=True)
        print("Done")

        print("Calculating hashes for tree 2...")
        hashes_tree2 = add_hashes(tree2)
        print("Done")

        print("Comparing trees...")
        if args.posterior_threshold:
            heights, subtrees = get_common_subtrees(tree1, subtree_times1, hashes_tree2,float(args.posterior_threshold))
        else:
            heights, subtrees = get_common_subtrees(tree1, subtree_times1, hashes_tree2)
        print("The median height of common subtrees is {}".format(round(np.median(heights),4)))
        
        with open(os.path.join(os.getcwd(),tree1_name + '_' + tree2_name + "_commontrees.txt"), 'w') as file:
            file.write("\n".join(subtrees) + "\n")
        file.close()

        with open(os.path.join(os.getcwd(), tree1_name + '_' + tree2_name + "_heights.txt"), 'w') as file:
            file.write("\n".join([str(h) for h in heights]) + "\n")
        file.close()


