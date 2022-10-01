# a function to list mean prevalence and total abundance on specified taxon level

function(psobject, taxalevel) {

  prevalencedf = apply(X = phyloseq::otu_table(psobject),
                       MARGIN = 1,
                       FUN = function(x){sum(x > 0)})   # how often  the ASV appears (sum of presence)

  prevalencedf = data.frame(Prevalence = prevalencedf,
                            TotalAbundance = taxa_sums(psobject),  # sum of abundance
                            as.data.frame(phyloseq::tax_table(psobject)))

  prevalencedf <-  prevalencedf %>%
    dplyr::group_by( {{taxalevel}} ) %>%
    dplyr::summarise(Mean_prevalence = mean(Prevalence),
                     Total_abundance = sum(TotalAbundance)) %>%
    dplyr::mutate(Rel_abundance = (Total_abundance / sum(Total_abundance) *100)) %>%
    dplyr::mutate(Cumulated = cumsum(Rel_abundance))

  count.taxon <- phyloseq::tax_table(psobject)@.Data %>%
    as.data.frame %>%
    rownames_to_column("ASVID") %>%
    dplyr::select(ASVID, {{taxalevel}}) %>%
    column_to_rownames("ASVID") %>%
    count({{taxalevel}})

  taxalevel <- enquo(taxalevel)
  prevalencedf <- prevalencedf %>%
    left_join(count.taxon, by = rlang::quo_name(taxalevel) ) %>%
    dplyr::arrange(desc(n))

return(prevalencedf)
}

