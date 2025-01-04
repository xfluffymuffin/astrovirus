library('randomcoloR')
library('colorspace')
library('dplyr')


#' Adds colors for taxa labels to metadata.
#' 
#' @param order_files  Path to text file with paths to files with taxa order (=order files).
#' If there are several text paths, they should be separated by \n.
#' The taxa in order files should be separated by \n.
#' 
#' @param metadata Path to text file with metadata for taxa (csv)
#' @return metadata dataframe with additional 'color' column 
add_colors2meta = function(order_files, metadata){
  
  meta = read.csv(metadata)

  # read file with paths to order files
  order_files = read.table(order_files)["V1"]
  list_taxa_df = list()
  for (i in 1:nrow(order_files)){
    taxa = read.table(order_files[i,1])
    colnames(taxa) = "GBAC"
    list_taxa_df = append(list_taxa_df,taxa)
  }
  
  # distinct colors (hex)
  colors = distinctColorPalette(length(list_taxa_df))
  print(colors)

  for (i in 1:length(list_taxa_df)){
    num_colors = length(list_taxa_df[[i]])
    # lightened color
    print(i)
    print(colors[i])
    color1 = lighten(colors[i], 0.4)
    # darkened color
    color2 = darken(colors[i], 0.4)
    
    print(color1)
    print(color2)
    # create pallete
    color_range <- colorRampPalette(c(color1, color2))
    # get range of colors
    colors_clade <- color_range(num_colors)
    
    # create dataframe of colors and corresponding names
    df = as.data.frame(cbind(list_taxa_df[[i]], colors_clade))
    colnames(df) = c("GBAC","color")
    
    # update list
    list_taxa_df[[i]] = df
  }
  # combine list of dataframes into one dataframe by row
  df_colors = bind_rows(list_taxa_df)
  meta_upd = full_join(meta, df_colors, by="GBAC")
  return(meta_upd)
}
