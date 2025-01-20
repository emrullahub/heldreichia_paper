# This script performs Gene Ontology (GO) enrichment analysis on a list of genes, focusing on 
# Biological Processes (BP). It also visualizes the results and generates detailed output files for 
# enriched pathways and gene-pathway associations. Key steps include:
# 1. Reading the gene list and performing GO enrichment analysis using the org.At.tair.db database.
# 2. Visualizing the results using barplot, dotplot, upset plot, and network plot.
# 3. Extracting and saving the enriched pathways and associated genes in a structured format.

library(tidyverse)
library(RColorBrewer)
library(clusterProfiler)
library(DOSE)
library(pheatmap)
library(enrichplot)
library(ggupset)
library(org.At.tair.db)
library(dplyr)
library(dplyr)

uniq_genes <- read.table("C:/Users/Pc/Documents/rcodes/randomforest_denemeler/geneontology/allgenes.only", stringsAsFactors = FALSE)$V1
go_enrichment <- enrichGO(gene = uniq_genes,
                          OrgDb = org.At.tair.db,
                          keyType = "TAIR",
                          ont = "BP",  # Ontology type: BP (Biological Process), CC (Cell Component), MF (Molecular Function)
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2)

barplot(go_enrichment, showCategory = 20, title = "GO Enrichment (BP)", font.size = 10)
dotplot(go_enrichment, showCategory = 20, title = "GO Enrichment (BP)")
upsetplot(go_enrichment, n = 20)

go_enrichment_plots <- pairwise_termsim(go_enrichment)
emapplot(go_enrichment_plots, showCategory = 20)
cnetplot(go_enrichment_plots, showCategory = 20)


go_results <- as.data.frame(go_enrichment)
#gene sets and pvalues
go_genes_pvalues <- go_results[, c("ID", "Description", "pvalue", "qvalue", "geneID")]
head(go_genes_pvalues)


gene_pathways <- go_enrichment@result %>%
  dplyr::select(geneID, Description) %>%
  tidyr::separate_rows(geneID, sep = "/") %>%
  dplyr::rename(Gene = geneID, Pathway = Description)
print(gene_pathways, n=Inf)
write.table(gene_pathways, file = "gene_pathways_table.txt", sep = "\t", row.names = FALSE, quote = FALSE)


gene_pathways_combined <- gene_pathways %>%
  group_by(Gene) %>%
  summarise(Pathways = paste(unique(Pathway), collapse = ", ")) %>%
  ungroup()
print(gene_pathways_combined, n = Inf)
write.table(gene_pathways_combined, file = "C:/Users/Pc/Desktop/heldbitis/gene_pathways_combined.txt", sep = "\t", row.names = FALSE, quote = FALSE)



enriched_genes <- unique(gene_pathways$Gene)
missing_genes <- setdiff(uniq_genes, enriched_genes)
print(missing_genes)