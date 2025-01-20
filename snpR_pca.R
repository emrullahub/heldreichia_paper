# This script performs PCA (Principal Component Analysis) on SNP data to visualize the population structure
# based on phenotypes and geographical distribution. Key steps include:
# 1. Reading and processing VCF data to extract SNP genotypes and metadata.
# 2. Converting genotypes into compatible formats (snpR object) for PCA and clustering analysis.
# 3. Performing PCA and extracting population centers and loadings for visualization.
# 4. Creating enhanced PCA scatter plots to depict population structure:

library(snpR)
library(vcfR)
library(ranger)
library(adegenet)
library(dplyr)
library(ggplot2)
library(ggrepel)

vcf_data <- read.vcfR("/home/emr/heldreichia/filtered_snpslast.vcf.gz")#
gll=vcfR2genlight(vcf_data)
x=as.matrix(gll)
gi=as.genind(x)
genotypes <- extract.gt(vcf_data, element = "GT")
alt_alleles <- vcf_data@fix[, "ALT"]
ref_alleles <- vcf_data@fix[, "REF"]
chr_info <- vcf_data@fix[, "CHROM"]
position_info <- vcf_data@fix[, "POS"]
convert_genotypes <- function(genotypes, ref_alleles, alt_alleles, chr_info, position_info) {
  genotype_matrix <- matrix(NA, nrow = nrow(genotypes), ncol = ncol(genotypes))
  
  for (i in 1:nrow(genotypes)) {
    for (j in 1:ncol(genotypes)) {
      gt <- genotypes[i, j]
      
      if (is.na(gt)) {
        genotype_matrix[i, j] <- "NN"
      } else if (gt == "0/0") {
        genotype_matrix[i, j] <- paste0(ref_alleles[i], ref_alleles[i])
      } else if (gt == "1/1") {
        genotype_matrix[i, j] <- paste0(alt_alleles[i], alt_alleles[i]) 
      } else if (gt == "0/1" || gt == "1/0") {
        genotype_matrix[i, j] <- paste0(ref_alleles[i], alt_alleles[i])
      } else {
        genotype_matrix[i, j] <- "NN"
      }
    }
  }
  genotype_df <- data.frame(Chr = chr_info, Position = position_info, genotype_matrix)
  colnames(genotype_df)[3:ncol(genotype_df)] <- colnames(genotypes)
  rownames(genotype_df) <- 1:nrow(genotype_df)
  
  return(genotype_df)
}
genotype_df <- convert_genotypes(genotypes, ref_alleles, alt_alleles, chr_info, position_info)
rownames(genotype_df) <- paste0(chr_info, "_", position_info)
genotypes <- as.data.frame(genotype_df)
genos <- genotypes[,-c(1:2)]
snp_meta <- genotypes[,1:2]
meta_data <- read.table("/home/emr/heldreichia/meta.list", header = FALSE)#
sample_meta <- data.frame(pop = meta_data[[2]], phenotype = meta_data[[3]], geography = meta_data[[4]], stringsAsFactors = FALSE)

snps <- import.snpR.data(genos, snp.meta = snp_meta, sample.meta = sample_meta, mDat = "NN")

getdata <- plot_clusters(snps, "pop.phenotype")
pca_data <- as.data.frame(getdata$data$pca)
pca_loadings <- as.numeric(getdata$pca_loadings)
pop_centers <- pca_data %>%
  group_by(pop) %>%
  summarize(
    PC1 = mean(PC1),
    PC2 = mean(PC2),
    phenotype = first(phenotype)  # Fenotip bilgisini ekle
  )





# type 1
ggplot(pca_data, aes(x = PC1, y = PC2, color = phenotype, shape = geography)) +
  geom_point(size = 6, alpha = 0.8, stroke = 1.2) +  # Şekillerin dış hatları siyah
  labs(
    title = NULL,
    x = paste0("PC1 (", round(pca_loadings[1], 2), "%)"),
    y = paste0("PC2 (", round(pca_loadings[2], 2), "%)")
  ) +
  theme_minimal(base_size = 18) +  # Temel yazı boyutu artırıldı
  theme(
    axis.title = element_text(size = 18, face = "bold"),  # Daha büyük ve kalın eksen başlıkları
    axis.text = element_text(size = 16),  # Daha büyük eksen değerleri
    legend.title = element_text(size = 18, face = "bold"),  # Daha büyük lejant başlığı
    legend.text = element_text(size = 16),  # Daha büyük lejant metni
    panel.grid.major = element_blank(),  # Ana grid çizgileri kaldırıldı
    panel.grid.minor = element_blank(),
    axis.line = element_line(linewidth = 1.2)  # Daha belirgin eksen çizgileri
  ) +
  scale_color_manual(
    values = c("a" = "#F4A582", "b" = "#92C5DE", "c" = "#B2182B")  # Renkler phenotype için
  ) +
  scale_shape_manual(
    values = c("x" = 16, "y" = 17, "z" = 18)  # Şekiller coğrafyaya göre
  )



# type 2
ggplot(pca_data, aes(x = PC1, y = PC2, fill = phenotype, shape = geography)) +
  geom_point(size = 6, alpha = 0.8, stroke = 1.2, color = "black") +  # Dış hatlar siyah, iç renkler fill ile
  labs(
    title = NULL,
    x = paste0("PC1 (", round(pca_loadings[1], 2), "%)"),
    y = paste0("PC2 (", round(pca_loadings[2], 2), "%)")
  ) +
  theme_minimal(base_size = 18) +  # Temel yazı boyutu artırıldı
  theme(
    axis.title = element_text(size = 18, face = "bold"),  # Daha büyük ve kalın eksen başlıkları
    axis.text = element_text(size = 16),  # Daha büyük eksen değerleri
    legend.title = element_text(size = 18, face = "bold"),  # Daha büyük lejant başlığı
    legend.text = element_text(size = 16),  # Daha büyük lejant metni
    panel.grid.major = element_blank(),  # Ana grid çizgileri kaldırıldı
    panel.grid.minor = element_blank(),
    axis.line = element_line(linewidth = 1.2)  # Daha belirgin eksen çizgileri
  ) +
  scale_fill_manual(
    values = c("a" = "#F4A582", "b" = "#92C5DE", "c" = "#B2182B")  # Phenotype renkleri
  ) +
  scale_shape_manual(
    values = c("x" = 21, "y" = 24, "z" = 22)  # Şekiller coğrafyaya göre
  )
