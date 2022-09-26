
pf.tree.with.legend <- function(pf, layout = "circular", factorvector = NULL) {

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
gg <- ggtree::ggtree(pf$tree,layout=layout, aes(alpha = 0.3))
gg <- gg+ggtree::geom_hilight(data = df, aes(node=node, fill=Phylofactor))

return(gg)
}
