##### Function for randomly selecting colors from two sets of colours.

mytreecolorfunction <- function(phyloseqobj, brightcols, greycols, prevelancethreshold, seed = 123) {      #orderfactor = the explanatory variable used in model
  # extract bins from phylofactor object
  prevelancedf = apply(X = phyloseq::otu_table(phyloseqobj),   
                       MARGIN = 1,
                       FUN = function(x){sum(x > 0)})
  prevelancedf = data.frame(Prevalence = prevelancedf,          
                            TotalAbundance = taxa_sums(phyloseqobj),
                            phyloseq::tax_table(phyloseqobj))
  prevelancedf <- plyr::ddply(prevelancedf, "Phylum", function(df1){
    data.frame(mean_prevalence=mean(df1$Prevalence),total_abundance=sum(df1$TotalAbundance,na.rm = T),stringsAsFactors = F)
  })
  
  # Create a vector with colors, depending on prevalence
    set.seed(seed)
    phyla <- nrow(prevelancedf)   # loop number
    colvec <- vector("character")   # open a vector
      for(phylum in 1:phyla)  {   # for loop to fill vector
        boom <- ifelse( ( prevelancedf[phylum,]$total_abundance > prevelancethreshold), base::sample(brightcols, size = 1), base::sample(greycols, size = 1))
        colvec <- c(colvec, boom)
  }
  names(colvec) <- prevelancedf$Phylum
  zerovec <- c("white")
  names(zerovec) <- "No colour"
  colvec <- c(zerovec, colvec)
  return(colvec)
}

