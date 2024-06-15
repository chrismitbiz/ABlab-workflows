function(ps.ILR,
         factorvector,
         explanatory = NULL,
         labelvec,
         ncol = 2) {

# plot ILRs
# Extract plot data and pivot longer by ILR factors
plotdata <- data.frame(sample_data(ps.ILR)) %>%
  pivot_longer(cols = starts_with("Phylofactor"), names_to = "ILR")
# Make ILRs leveled factors so they appear in right order
plotdata$ILR <- factor(plotdata$ILR, levels = plotdata[1:PF$nfactors,]$ILR)
# Create colours (same as the pylofactors above) and name color vector by ILR name
# the name of the vector are required so ggpolot automatically allocated each colour to the right ILR
# cols <- viridis::turbo(nrow(PF$factors), direction = 1)
# names(cols) <- plotdata[1:PF$nfactors,]$ILR

plotdata_sub <- plotdata %>%
  filter(ILR %in% factorvector)
# ggplot with ggpubr
plotdata_sub$ILR <- factor(plotdata_sub$ILR,
                           levels = plotdata_sub[1:length(unique(plotdata_sub$ILR)),]$ILR)

plotdata_sub %>%
  ggboxplot(x = explanatory, y = "value",
            combine = TRUE, facet.by = "ILR",
            fill = "ILR", ncol = ncol, alpha = 0.5,
            ylab = ("Isometric log ratio (ILR)"),
            panel.labs.font = list(size = 11),
            panel.labs = list(ILR = labelvec) )  +
  geom_point(position = "jitter", alpha = 0.5) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12) ) +
  xlab(bquote(CO[2]) )


}
