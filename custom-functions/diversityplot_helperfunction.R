# Helper function to plot diversity indices into a given phyloseq object
# Default = Observed, change to any column in sample_data(psobject)

diversityplot_helperfunction <- function(psobject, y = "Observed") {
  
  p <- data.frame(phyloseq::sample_data(psobject)) %>%
    ggpubr::ggboxplot(., x = "Digester", 
                      y = y, 
                      fill = "Digester", 
                      ylab = paste(y), 
                      xlab = "Digester") +
    theme_bw() +
    theme( panel.grid.major = element_blank(), 
           panel.grid.minor = element_blank() ) +
    geom_jitter(width = 0.1) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    #    theme(axis.title.x=element_blank(),  # Options
    #        axis.text.x=element_blank(),
    #        axis.ticks.x=element_blank()) +
    ggpubr::stat_compare_means(method = "kruskal",label.x = 2) 
  p
}
