#Finds taxa names in string in newick format
find_taxa = function(stri){

  matches = gregexpr("[A-Za-z0-9_./]+",stri)
  
  if (matches[[1]][1] != -1 ){
    return(regmatches(stri, matches)[[1]])
  }
  else{
    return(c())
  }
}

# Produces list of taxa names in each coinciding subtree

get_subtrees = function(file_subtrees){
  lines = readLines(file_subtrees)
  names_taxa = lapply(lines,find_taxa)
  names(names_taxa) <- rep("same", length(names_taxa))
  return(names_taxa)
}
