# This script performs a random forest analysis on SNP data to identify significant SNPs 
# associated with a phenotype, using permutation tests for significance thresholds. 
# Key steps include:
# 1. Reading and processing VCF data to extract genotypes and SNP metadata.
# 2. Converting genotype data into formats compatible with SNP and random forest analysis.
# 3. Performing random forest classification to evaluate SNP importance for phenotype prediction.
# 4. Running permutation tests to establish significance thresholds for SNP importance.
# 5. Extracting significant SNPs based on Mean Decrease Accuracy (MDA) and generating output files.
# 6. Matching significant SNPs to their chromosome and position for identifying associated genes, 
#    which can then be used as input for heatmap visualization in scripts like dapcnet.R.

library(snpR)
library(vcfR)
library(ranger)
library(adegenet)
library(dplyr)
library(ggplot2)
library(genepop)

vcf_data <- read.vcfR("C:/Users/Pc/Desktop/snps_filtered.vcf.gz")
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
meta_data <- read.table("C:/Users/Pc/Downloads/meta.list", header = FALSE)#
sample_meta <- data.frame(
  pop = meta_data[[2]],     # 2. column: population
  phenotype = meta_data[[3]],       # 3. column: phenotype
  stringsAsFactors = FALSE
)
snps <- import.snpR.data(genos, snp.meta = snp_meta, sample.meta = sample_meta, mDat = "NN")
is.snpRdata(snps)


set.seed(123) # one run
rf <- run_random_forest(snps, response = "phenotype", pvals = FALSE, num.trees = 10000)


#func
run_permutation_test <- function(snp_data, response_col, nperm = 1000, num_trees = 10000) {
  # Orijinal model
  set.seed(123)
  original_rf <- run_random_forest(snp_data, response = response_col, pvals = FALSE, num.trees = num_trees)
  original_importance <- original_rf$models$.base_.base$model$variable.importance
  
  # Genotip matrisini al
  genotype_matrix <- snp_data@geno.tables$gs
  
  # Permütasyon için boş matris
  null_importance <- matrix(NA, nrow = length(original_importance), ncol = nperm)
  rownames(null_importance) <- names(original_importance)
  
  # Permütasyon döngüsü
  for (i in 1:nperm) {
    # Genotip verilerini karıştır
    permuted_genotype_matrix <- apply(genotype_matrix, 2, sample)
    
    # Karıştırılmış genotip verileriyle yeni bir `snpRdata` nesnesi oluştur
    permuted_snp_data <- snp_data
    permuted_snp_data@geno.tables$gs <- permuted_genotype_matrix
    
    # Karıştırılmış SNP verileriyle rastgele f çalıştır
    permuted_rf <- run_random_forest(permuted_snp_data, response = response_col, pvals = FALSE, num.trees = num_trees)
    
    # Permütasyon önem derecelerini kaydet
    null_importance[, i] <- permuted_rf$models$.base_.base$model$variable.importance
  }
  
  # Sonuçları döndür
  return(list(original_importance = original_importance, null_importance = null_importance))
}
permutation_results <- run_permutation_test(snps, "phenotype", nperm = 1000, num_trees = 10000)


null_importance_means <- rowMeans(permutation_results$null_importance, na.rm = TRUE)
null_importance_threshold <- quantile(null_importance_means, 0.99)

significant_snps <- names(permutation_results$original_importance[
  permutation_results$original_importance > null_importance_threshold
])
print(paste("important snp count:", length(significant_snps)))
print(null_importance_threshold)



null_importance_mda <- apply(permutation_results$null_importance, 1, mean, na.rm = TRUE)
write.table(null_importance_mda, file = "null_importance_mda_table.txt", 
            row.names = TRUE, col.names = TRUE, quote = FALSE)
mda_threshold <- quantile(null_importance_mda, 0.99)


original_importance <- permutation_results$original_importance
write.table(original_importance, file = "C:/Users/Pc/Desktop/outputs/original_importance_table.txt", 
            row.names = TRUE, col.names = TRUE, quote = FALSE)
mda_values <- original_importance

mda_threshold <- quantile(null_importance_means, 0.99)

significant_snps <- names(original_importance[mda_values > mda_threshold])
significant_snps_mda <- original_importance[significant_snps]
write.table(significant_snps_mda, file = "C:/Users/Pc/Desktop/outputs/significant_snps_mda_table.txt", 
            row.names = TRUE, col.names = TRUE, quote = FALSE)
result_table <- data.frame(SNP = significant_snps, MDA = significant_snps_mda)

print(result_table)

original_importance[mda_values > mda_threshold] # print SNPs and MDA 

write.table(result_table, file = "result_table.txt", 
            row.names = TRUE, col.names = TRUE, quote = FALSE)



# sort it!
sorted_snps <- sort(significant_snps, decreasing = TRUE)
print(sorted_snps)



sorted_result_table <- result_table[order(-result_table$MDA), ]

print(sorted_result_table)
write.table(sorted_result_table, file = "sorted_result_table.txt", 
            row.names = TRUE, col.names = TRUE, quote = FALSE)



hist(result_table$MDA,
     breaks = 10,              
     col = "blue",              
     border = "black",          
     main = "MDA ",
     xlab = "MDA")


vcf_snp_positions <- data.frame(
  CHROM = vcf_data@fix[, "CHROM"],
  POS = as.numeric(vcf_data@fix[, "POS"])
)
matched_positions <- vcf_snp_positions[vcf_snp_positions$POS %in% sorted_snps, ]
head(matched_positions)
write.table(matched_positions, file = "C:/Users/Pc/Desktop/outputs/matched_snp_positions.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE)





sorted_result_table <- read.table("sorted_result_table_04.txt", header = TRUE)
sorted_snps <- sorted_result_table$SNP
vcf_snp_positions <- data.frame(
  CHROM = vcf_data@fix[, "CHROM"],
  POS = as.numeric(vcf_data@fix[, "POS"])
)
matched_positions <- vcf_snp_positions[sorted_snps, ]
matched_positions$SNP_satiri <- sorted_snps
matched_positions <- matched_positions[, c("SNP_satiri", "CHROM", "POS")]
write.table(matched_positions, file = "C:/Users/Pc/Desktop/outputs/SNP_positions_output.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
