## addILRtophylose function ##
# Function to take a phylofactor file and extract the mean sample ILRs and add those to the phyloseq sample data
# those can then be used for plotting graphs together with abundance data.
# Add name of phyloseqobject and phylofactorobject

function(phyloseqobject, phylofactorobject) {
  # Export taxa df from phyloseq object
  taxadf <- as.data.frame(tax_table(phyloseqobject)@.Data)
  taxadf <- unite(taxadf, taxonomy, sep = ";", remove = TRUE, na.rm = FALSE) %>% rownames_to_column("OTU_ID")
  # export metadata from phyloseq object
  meta.df <- data.frame(sample_data(phyloseqobject)) %>% rownames_to_column("OTU_ID")
  # add ILRs for each sample and factor to the metadataframe
  factornumber <- nrow(phylofactorobject$factors)
  for (i in 1:factornumber) {
    smry <- phylofactor::pf.summary(phylofactorobject, taxadf,  factor=i)
    ilr.df <- data.frame(mean.ilr = smry$ilr) %>% rownames_to_column("OTU_ID")
    meta.df <- meta.df %>% left_join(ilr.df, by = "OTU_ID")
    names(meta.df)[ncol(meta.df)] <- paste0("Phylofactor_", i)
  }
  # re-import metadata into phyloseq object
  phyloseqobject@sam_data <- phyloseq::sample_data(meta.df)
  return(phyloseqobject)
}
