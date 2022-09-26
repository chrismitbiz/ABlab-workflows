## phyfacsummary function ##
## Function to summaries all phylofactor factors in a meaningful way into a dataframe
## Enter the phylofactor file and the phyloseq object to get a dataframe with model results, abundances and taxa info
## the third function input is '4' or '5',  4 is used when the phylofactor model setting was set to 'F', if 'var' use 5
## p adjusted method = Bonferroni

phyfacsummary <- function(phylofactorfile, phyloseqobject, Var_or_F) {
  datalist = list()
  factornumber <- nrow(phylofactorfile$factors)
  taxadf <- as.data.frame(tax_table(phyloseqobject)@.Data)
  taxadf <- unite(taxadf, taxonomy, sep = ";", remove = TRUE, na.rm = FALSE) %>% rownames_to_column("OTU")
  otudf <- data.frame(Abundance = rowSums(as.data.frame(otu_table(phyloseqobject)@.Data))) %>% rownames_to_column("OTU")
  for (i in 1:factornumber) {
    smry <- phylofactor::pf.summary(phylofactorfile, taxadf,  factor=i)
    td <- phylofactor::pf.tidy(smry)
    factorsgroup <- data.frame(smry$group1) %>% dplyr::select(IDs.otuIDs,IDs.TaxaIDs) %>% 
      dplyr::rename(OTU = IDs.otuIDs) %>% dplyr::rename(Taxa_ID = IDs.TaxaIDs) %>% 
      mutate(Phylofactor = i) %>% 
      mutate('Pr(>F)' = phylofactorfile$factors[[i,Var_or_F]] )
      ncoefs <- nrow(data.frame(td$Coefficients))
        for(c in 1:ncoefs) {                                   # Head of for-loop
          new <- rep(td$Coefficients[[c]], nrow(factorsgroup)) # Create new column
          factorsgroup[ , ncol(factorsgroup) + 1] <- new      # Append new column
          colnames(factorsgroup)[ncol(factorsgroup)] <- paste0(rownames(data.frame(td$Coefficients[c])))  # Rename column name
          }
  datalist[[i]] <- factorsgroup
  }
  
phylofactorsummary <- dplyr::bind_rows(datalist) 
  phylofactorsummary$`Pr(>F)adj.` <- p.adjust(phylofactorsummary$`Pr(>F)`, "fdr")
  return(phylofactorsummary %>% left_join(otudf, by = "OTU") %>% as_tibble() %>% separate(Taxa_ID, sep=";", c("Kingdom","Phylum","Class","Order","Family","Genus","Species")))
}
