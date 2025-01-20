# This script performs DAPC (Discriminant Analysis of Principal Components) analysis on SNP data extracted 
# from a VCF file. The code processes the data to identify population structure and visualize SNP 
# contributions and membership probabilities. Key steps include:
# 1. Reading and converting VCF data into genind and genlight formats for DAPC.
# 2. Associating phenotypic information (population groups) with samples.
# 3. Performing DAPC and visualizing results using scatter plots and barplots.
# 4. Using SNPzip to identify key SNPs contributing to discrimination between groups.
# 5. Generating a heatmap to display genotypes of selected SNPs across individuals.

library(vcfR)
library(adegenet)
library(ade4)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

vcf <- read.vcfR("filtered_snps.vcf.gz")
gll=vcfR2genlight(vcf)
x=as.matrix(gll)
gi=as.genind(x)

sample_list <- read.table("sample.list", header = FALSE, stringsAsFactors = FALSE)
# list sample
#bou_alakir_5658_01 a bou_alakir
#bou_gokbel_5673_01 a bou_gokbel
#bou_makdag_5643_01 a bou_makdag
#rot_erentp_5840_08 c rot_erentp
individual_names <- sample_list$V1
phenotype <- sample_list$V2
r <- factor(phenotype)


dp=dapc(gi, var.contrib = T, pop=r, perc.pca=80, scale = FALSE, n.da = 1)


my_colors <- c("#F4A582", "#92C5DE", "#B2182B")
scatter(dp, col = my_colors, pch = 19, scree.da = TRUE, legend = TRUE)


compoplot(dp, col = my_colors)
# or
assignments <- dp$posterior
par(mar = c(6, 4, 10, 8), xpd = TRUE)# margin
bp <- barplot(t(assignments), beside = FALSE, col = my_colors, border = NA, space = 0.4, xaxt = "n", yaxt = "n")
mtext(text = individual_names, side = 3, line = 1, las = 2, cex = 0.8, col = "black", at = bp - 0.3)
axis(2, las = 1, cex.axis = 0.8, col = "black", col.axis = "black")
mtext("Membership Probability", side = 2, line = 2, cex = 1.1, col = "black")
mtext("Individuals", side = 1, line = 0.4, cex = 1.1, col = "black")
legend("bottomright", inset = c(0, -0.1), legend = c("a", "b", "c"), fill = my_colors, 
       horiz = TRUE, cex = 1.1, bty = "n", x.intersp = 0.5, y.intersp = 1.2)
box(lwd = 2)


snpzip_results <- snpzip(gi, dp, plot = FALSE, xval.plot = FALSE, loading.plot = TRUE, method = "average")

membership_matrix <- dp$posterior
assigned_groups <- apply(membership_matrix, 1, which.max)
group_names <- colnames(membership_matrix)[assigned_groups]
group_counts <- table(group_names)
print(group_counts)




#selected_snps <- read.table("randomforest_snps.txt")$V1 # random forest data

positions <- snpzip_results$FS$`Contributions of selected alleles to discriminant axis`
selected_snps <- names(positions) # dapc data

pop_info <- sample_list$V3
row_annotation <- data.frame(Population = pop_info)
rownames(row_annotation) <- indNames(gi)

genotype_matrix <- tab(gi)
selected_snp_matrix <- genotype_matrix[, selected_snps, drop = FALSE]

row_names <- indNames(gi)  
col_names <- selected_snps   

heatmap_colors <- colorRampPalette(brewer.pal(9, "YlOrRd"))(12)
pheatmap(
  selected_snp_matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = heatmap_colors,
  main = "Selected SNPs Heatmap",
  fontsize = 8,
  fontsize_row = 5,
  fontsize_col = 7,
  show_colnames = TRUE,
  show_rownames = TRUE,
  labels_row = row_names,
  labels_col = col_names,
  annotation_row = row_annotation
)