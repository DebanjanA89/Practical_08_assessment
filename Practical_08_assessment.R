#Installing and loading required packages
install.packages("BiocManager")
library(BiocManager)
install.packages('tidyverse')
library(tidyverse)
install(c('tximport', 'DESeq2', 'biomaRt', 'pheatmap'))
library(tximport)
library(DESeq2)
library(biomaRt)
library(pheatmap)
library(tidyverse)

sample_table = read_csv('https://raw.githubusercontent.com/sjcockell/mmb8052/main/practicals/practical_08/data/sample_table.csv')
files = pull(sample_table, Run)
files = paste0('counts/', files, '/quant.sf')
names(files) = pull(sample_table, Run)
gene_map = read_csv('https://github.com/sjcockell/mmb8052/raw/main/practicals/practical_08/extdata/gene_map.csv')
txi = tximport(files,
               type='salmon',
               tx2gene=gene_map,
               ignoreTxVersion=TRUE)

getwd()
#To set the correct working directory
setwd("C:/Users/c2052841/OneDrive - Newcastle University/debanjana")
txi = tximport(files,
               type='salmon',
               tx2gene=gene_map,
               ignoreTxVersion=TRUE)
gene_map = read_csv('https://github.com/sjcockell/mmb8052/raw/main/practicals/practical_08/extdata/gene_map.csv')
dds = DESeqDataSetFromTximport(txi, colData = sample_table, design = ~ Group)
dds = estimateSizeFactors(dds)
dds = estimateDispersions(dds)
dds = estimateDispersions(dds)
dds = nbinomWaldTest(dds)
plotDispEsts(dds)

#PCA plot and heatmap
rld = rlog(dds)
plotPCA(rld, intgroup='Group')
library(ggplot2)
plotPCA(rld, intgroup='Group') + theme_bw( )
sample_distance = dist(t(assay(rld)), method='euclidian')
sample_distance_matrix = as.matrix(sample_distance)
> heatmap_annotation = data.frame(group=colData(dds)[,c('Group')], row.names=rownames(colData(dds)))
heatmap_annotation = data.frame(group=colData(dds)[,c('Group')], row.names=rownames(colData(dds)))

#Installing RColorBrewer to improve heatmap
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


results_table = results(dds, contrast= c('Group', 'Allo24h', 'Naive'))
summary(results_table)
library(dplyr)
results_tibble = as_tibble(results_table, rownames='ensembl_gene_id')
filtered_results = filter(results_tibble, complete.cases(results_tibble))
filtered_results = mutate(filtered_results, logPVal = -log10(padj))
ggplot(filtered_results, aes(x=log2FoldChange, y=logPVal)) +
  +   geom_point(aes(colour=logPVal)) +
  +   geom_vline(xintercept=1) +
  +   geom_vline(xintercept= -1) +
  +   geom_hline(yintercept= -log10(0.05), linetype=2, colour="dodgerblue")

#Repeating the same code for Allo2h
results_table1 = results(dds, contrast= c('Group', 'Allo2h', 'Naive'))
summary(results_table1)
results_tibble1 = as_tibble(results_table1, rownames='ensembl_gene_id')
filtered_results1 = filter(results_tibble1, complete.cases(results_tibble1))
filtered_results1 = mutate(filtered_results1, logPVal = -log10(padj))
ggplot(filtered_results1, aes(x=log2FoldChange, y=logPVal)) +
  +   geom_point(aes(colour=logPVal)) +
  +   geom_vline(xintercept=1) +
  +   geom_vline(xintercept= -1) +
  +   geom_hline(yintercept= -log10(0.05), linetype=2, colour="dodgerblue")

ensembl108 = useEnsembl(biomart="ensembl", version=108)
library(biomaRt)
ensembl108 = useEnsembl(biomart="ensembl", version=108)
ensembl108 = useDataset("mmusculus_gene_ensembl", mart=ensembl108)
annotation = getBM(attributes=c('ensembl_gene_id', 'chromosome_name',
                                'start_position', 'end_position',
                                'strand', 'gene_biotype', 'external_gene_name',
                                'description'),
                   filters = 'ensembl_gene_id', values = filtered_results$ensembl_gene_id,
                   mart = ensembl108)
annot_results = left_join(filtered_results, annotation)
annot_results = arrange(annot_results, padj)
View(head(annot_results, 10))
View(head(annot_results, 20))
degs = filter(annot_results, abs(log2FoldChange) > 1 & padj < 0.05)
View(head(degs, 10))
