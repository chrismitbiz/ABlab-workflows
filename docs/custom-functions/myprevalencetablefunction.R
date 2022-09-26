# a function to list mean prevalence and total abundance on specified taxon level

function(psobject, taxalevel) {

  prevalencedf = apply(X = phyloseq::otu_table(psobject),
                       MARGIN = 1,
                       FUN = function(x){sum(x > 0)})

  prevalencedf = data.frame(Prevalence = prevalencedf,
                            TotalAbundance = taxa_sums(psobject),
                            as.data.frame(phyloseq::tax_table(psobject)))

  prevalencedf <-  prevalencedf %>%
    dplyr::group_by( {{taxalevel}} ) %>%
    dplyr::summarise(Mean_prevalence = mean(Prevalence),
                     Total_abundance = sum(TotalAbundance)) %>%
    dplyr::mutate(Rel_abundance = (Total_abundance / sum(Total_abundance) *100)) %>%
    arrange(desc(Total_abundance)) %>%
    dplyr::mutate(Cumulated = cumsum(Rel_abundance))

return(prevalencedf)
}
