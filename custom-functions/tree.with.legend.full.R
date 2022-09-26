
pf.tree.with.legend.full <- function(pf, phyloseqfile, layout = "circular", factorvector = NULL) {

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

### add list grouping to tree and call it Phylum
df <- data.frame(factor.map, node = nd)
df$Phylofactor <- factor(df$Phylofactor, levels = df$Phylofactor)

# prepare for plotting
# Prevalence table to check which taxa to take out of the visualisation
prevelancedf = apply(X = phyloseq::otu_table(phyloseqfile),
                     MARGIN = 1,
                     FUN = function(x){sum(x > 0)})
prevelancedf = data.frame(Prevalence = prevelancedf,
                          TotalAbundance = taxa_sums(phyloseqfile),
                          phyloseq::tax_table(phyloseqfile))

prevelancedf <- plyr::ddply(prevelancedf, "Phylum", function(df1){
  data.frame(mean_prevalence=mean(df1$Prevalence), total_abundance=sum(df1$TotalAbundance,na.rm = T), stringsAsFactors = F)
})
# selection of phylum names to keep
legendlimit <- (prevelancedf %>% dplyr::filter(total_abundance > 100))$Phylum

# selection to filter out
losers <- (prevelancedf %>% dplyr::filter(total_abundance < 100))$Phylum

## Extract a taxadataframe based on your final filters
taxa.df <- phyloseq::tax_table(phyloseqfile)@.Data %>% as.data.frame
## Create a list to split otus by phyla, this will be used for colour-grouping the tree
taxa.df.phylum <- taxa.df %>%
  rownames_to_column("OTUID") %>%
  dplyr::select(OTUID, Phylum) %>%
  column_to_rownames("OTUID")
phyla.list <- split(rownames(taxa.df.phylum), taxa.df.phylum$Phylum, drop = TRUE)
## adding that phyla list to the tree
tree.df <-  phyloseq::phy_tree(pf$tree)
### add list grouping to tree and call it Phylum
ggtree_gps <- ggtree::groupOTU(tree.df, phyla.list, "Phylum")
### Plot tree
gg <- ggtree::ggtree(ggtree_gps,
               aes(color=Phylum, alpha = 0.3),
               layout = layout) +
  #xlim(-2.5, NA) +   # to shorten the tree tips
  theme(legend.position="right") +
  scale_colour_manual(values = treecols,
                      limits = legendlimit) +
  guides(color = guide_legend(override.aes = list(size=5)))

gg <- gg+ggtree::geom_hilight(data = df, aes(node=node, fill=Phylofactor))

return(gg)
}
