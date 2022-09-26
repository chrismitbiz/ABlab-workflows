
function(pf, phyloseqfile,
          layout = "circular",
          branch.length= "branch.length",
          factorvector = NULL,
          ordervector = NULL)     {

# Step one: Create a dataframe with phylofactor nodes to highlight
  factor.map=data.frame('Phylofactor'=1:pf$nfactors,'group'=rep(1,pf$nfactors))

  if (is.null(factorvector)) {
    factor.map=data.frame('Phylofactor'=1:pf$nfactors,'group'=rep(1,pf$nfactors))
  }  else {
    factor.map <- factor.map %>% dplyr::filter(Phylofactor %in% factorvector)
  }
  m <- nrow(factor.map)
  #color.fcn=viridis::viridis
  #cols <- color.fcn(m)
  #factor.map$colors <- cols

  n=ape::Ntip(pf$tree)
  method = "factors"
  nd <- NULL
  for (i in 1:m){
    if (method=='factors'){
      grp <- pf$tree$tip.label[pf$groups[[factor.map[i,1]]][[factor.map[i,2]]]]
    } else {
      grp <- pf$tree$tip.label[GroupList[[i]]]
    }

    grp <- intersect(grp, pf$tree$tip.label)
    if (length(grp)>0 & any(setdiff(pf$tree$tip.label,grp) %in% pf$tree$tip.label)){
      if (length(grp)>1){
        nd <- c(nd, tidytree::MRCA(pf$tree,grp))
      } else {
        nd <- c(nd,match(grp,pf$tree$tip.label))
      }
    }
  }

## create dataframe with nodes to that identify the phylofactors
  df <- data.frame(factor.map, node = nd)
  df$Phylofactor <- factor(df$Phylofactor, levels = df$Phylofactor)


# Step 2: Add groups to ASVs into the tree that help associate ASVs with their taxon
  ## Create a list to split otus by phyla, this will be used for colour-grouping the tree
  ## Extract a taxadataframe
  taxa.df.phylum <- phyloseq::tax_table(phyloseqfile)@.Data %>%
    as.data.frame %>%
    rownames_to_column("OTUID") %>%
    dplyr::select(OTUID, Phylum) %>%
    column_to_rownames("OTUID")
 phyla.list <- split(rownames(taxa.df.phylum),
                     taxa.df.phylum$Phylum,
                     drop = TRUE)

  ## adding that phyla list to the tree
  tree.df <-  phyloseq::phy_tree(pf$tree)
  ### add list grouping to tree and call it Phylum
  ggtree_gps <- ggtree::groupOTU(tree.df,
                                 phyla.list, "Phylum")

# # Step 3: Plot tree
  ## order of phyla in legend and selection of phylum names to keep
  if (is.null(ordervector)) {
    ordervec <- (taxa.df.phylum %>%
                   count(Phylum) %>%
                   arrange(desc(n)) %>%
                   dplyr::filter(n > 2))$Phylum
  }  else {
    ordervec <- ordervector
  }


  gg <- ggtree::ggtree(ggtree_gps,
                       aes(color=Phylum, alpha = 0.3),
                       layout = layout,
                       branch.length = "none") +
    #xlim(-2.5, NA) +   # to shorten the tree tips
    theme(legend.position="right") +
    scale_colour_manual(values = treecols,
                        limits = ordervec,
                        breaks = ordervec) +
    guides(color = guide_legend(override.aes = list(size=5))) +
    guides(alpha = FALSE)

  gg <- gg+ggtree::geom_hilight(data = df, aes(node=node, fill=Phylofactor))

  return(gg)
}
