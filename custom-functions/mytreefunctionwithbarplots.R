# create a tree from a phylofactor object and a phyloseq object
# point is to highlight the phylofactors, show a colored tree and abundances for each tip
# need 'treecols' vector.

function(pf, ps.phylo,
         layout = "circular",
         branch.length = "none",
         factorvector = NULL,
         offset.text = 5,
         pwidth = 0.05) {

  # replace the tree with the unrooted tree from the phylofactor object
  phyloseq::phy_tree(ps.phylo) <- pf$tree
  # get overview of abundances, mean prevalence is the mean 'appearance' of ASVs of the taxon across all samples
  # knitr::kable(prevalencedf(ps.phylo, Phylum))

  melt_simple <- psmelt(ps.phylo) %>%
    #filter(Abundance  120) %>%
    dplyr::select(OTU, val=Abundance) %>% unique() #unquote metadatacolumn

  # Plot the tree
  ## First Add groups to ASVs into the tree that help associate ASVs with their taxon
  ## Create a list to split otus by phyla, this will be used for colour-grouping the tree
  ## Extract a taxadataframe
  taxa.df.phylum <- phyloseq::tax_table(ps.phylo)@.Data %>%
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

  p <-  ggtree::ggtree(ggtree_gps,
                       aes(color=Phylum),
                       layout = layout,
                       branch.length = "none") +
    #xlim(-2.5, NA) +   # to shorten the tree tips
    theme(legend.position="right") +
    scale_colour_manual(values = treecols,
                        limits = ordervec,
                        breaks = ordervec,
                        guide=guide_legend(override.aes = list(size=5),
                                           keywidth = 0.2,
                                           keyheight = 0.2,
                                           order=1))

  # Add barplots of abundances at the end
  #p <- rotate_tree(p, -90)
  p <-  p + ggnewscale::new_scale_fill() +
    geom_fruit(
      data=melt_simple,
      geom=geom_bar,
      mapping = aes(
        y=OTU,
        x=val,
        group=label,
        fill=Phylum),
      pwidth=0.38,
      orientation="y",
      stat="identity",
      show.legend = FALSE) +
    scale_fill_manual(values = treecols,
                      limits = ordervec,
                      breaks = ordervec)

  # Extract the phylofactor nodes (code copied from phylofactor package)
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
  ## add nodes highlights
  p <-  p +
    ggtree::geom_hilight(data = df,
                         aes(node=node),
                         fill="black",
                         alpha=0.3 )

  # add a clade annotation - identifying the phylofactor
  for(i in factorvector) {
    p <- p +  geom_cladelabel(node= (df %>% dplyr::filter(Phylofactor == i))$node,
                              label= paste("PF",(df %>% dplyr::filter(Phylofactor == i))$Phylofactor),
                              offset.text= offset.text,
                              hjust='center')
  }
  return(p)
}

