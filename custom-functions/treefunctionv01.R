# If only the tree is required, without highlights
treefunction <- function(phyloseqfile, layout = "circular") {
  ## Extract a taxadataframe based on your final filters
  taxa.df <- phyloseq::tax_table(phyloseqfile)@.Data %>% as.data.frame
  ## Create a list to split otus by phyla, this will be used for colour-grouping the tree
  taxa.df.phylum <- taxa.df %>%
    rownames_to_column("OTUID") %>%
    dplyr::select(OTUID, Phylum) %>%
    column_to_rownames("OTUID")
  phyla.list <- split(rownames(taxa.df.phylum), taxa.df.phylum$Phylum, drop = TRUE)
  ## adding that phyla list to the tree
  tree.df <-  phyloseq::phy_tree(phyloseqfile)
  ### add list grouping to tree and call it Phylum
  ggtree_gps <- ggtree::groupOTU(tree.df, phyla.list, "Phylum")
  ### Plot tree
  ggtree::ggtree(ggtree_gps,
                 aes(color=Phylum, alpha = 0.3),
                 layout = layout) +
    #xlim(-2.5, NA) +   # to shorten the tree tips
    theme(legend.position="right") +
    scale_colour_manual(values = treecols) +
    guides(color = guide_legend(override.aes = list(size=5)))
}
