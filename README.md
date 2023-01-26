# Practical_08_assessment

These are the codes that I have used for my project
I have used R studio version 4.2.1

I have made a few changes to the codes that were initially provided

# Change in the PCA Plot
plotPCA(rld, intgroup='Group') + theme_bw( ) 

#Change in heatmap color sequence and changed the colors allocated for Allo24h and 2h
install.packages("RColorBrewer")
display.brewer.all()
library(RColorBrewer)
display.brewer.all()
Purples = colorRampPalette(brewer.pal(9, "Purples"))(100)
library(pheatmap)
pheatmap(sample_distance_matrix,
         clustering_distance_rows=sample_distance,
         clustering_distance_cols=sample_distance,
         annotation_col = heatmap_annotation,
         color = Purples,
         annotation_colors = list ( 'group' = c(Allo24h ="#FF9999", Allo2h = "#00CC66", Naive = "#66B2FF")))
         
# Change in Volcano Plot
ggplot(filtered_results, aes(x=log2FoldChange, y=logPVal)) +
  +   geom_point(aes(colour=logPVal)) +
  +   geom_vline(xintercept=1) +
  +   geom_vline(xintercept= -1) +
  +   geom_hline(yintercept= -log10(0.05), linetype=2, colour="dodgerblue")
