# This script performs a random forest analysis on SNP data to classify phenotypes.
# The main steps are as follows:
# 1. Reads a CSV file containing genotype matrix and separates individual IDs, class labels, and SNP data.
# 2. Builds a random forest model using all SNPs and computes variable importance scores (impurity-based).
# 3. Filters SNPs above the 50th percentile of importance scores, rebuilds the random forest with the filtered SNPs, and saves updated importance scores.
# 4. For a range of importance percentiles, creates SNP subsets and evaluates each using random forest to calculate out-of-bag (OOB) error.
# 5. Identifies the SNP subset with the lowest OOB error and applies "backward purging": iteratively removes the least important SNP while tracking OOB error.
# 6. Saves the purged SNP list that yields the best OOB score. These SNPs can later be matched to genes and used in visualization tools such as heatmaps. (dapcnet.R)

# This code snippet is based on the original code developed by Çisel Kemahlı: https://github.com/ciselkemahli
# I would like to thank her for her contribution.

library(ranger)
csv_file <- "C:/Users/Pc/Desktop/randomf/ranger/class_data_noNA.csv"
class_data <- read.csv(csv_file, header = TRUE, stringsAsFactors = FALSE)

# Indv ID and Type
individuals <- class_data[, 1]   # indv ID
y <- as.factor(class_data[, 2])  # Type data
geno <- class_data[, -c(1, 2)]   # snps

# Importance
ntree <- 5000
mtry <- floor(sqrt(ncol(geno)))
set.seed(42)
cat("[1/5] Calculating importance with the full model...\n")
rf_full <- ranger(x = geno, y = y,
                  num.trees = ntree, mtry = mtry,
                  importance = "impurity", oob.error = TRUE)
imp_full <- sort(rf_full$variable.importance, decreasing = TRUE)
write.table(data.frame(SNP = names(imp_full), Importance = imp_full),
            file = "C:/Users/Pc/Desktop/randomf/ranger/OUTPUT/importance_full.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# %50 
thr <- quantile(imp_full, 0.5)
snps_prefilt <- names(imp_full)[imp_full >= thr]
write.table(snps_prefilt, file = "C:/Users/Pc/Desktop/randomf/ranger/OUTPUT/snps_prefilt.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# importance calc w filtered snps
rf_prefilt <- ranger(x = geno[, snps_prefilt, drop = FALSE], y = y,
                     num.trees = ntree, mtry = floor(sqrt(length(snps_prefilt))),
                     importance = "impurity", oob.error = TRUE)
imp_prefilt <- sort(rf_prefilt$variable.importance, decreasing = TRUE)
write.table(data.frame(SNP = names(imp_prefilt), Importance = imp_prefilt),
            file = "C:/Users/Pc/Desktop/randomf/ranger/OUTPUT/importance_prefilt.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# percentiles
percentiles <- c(0.98, 0.97, 0.96, 0.95, 0.90, 0.80, 0.70)
subset_results <- data.frame(Dilim = percentiles, NLokus = NA, OOB = NA)
write.table(subset_results, file = "C:/Users/Pc/Desktop/randomf/ranger/OUTPUT/subset_results.tsv", sep = "\t", row.names = FALSE)

cat("[3/6] Subset testing is starting...\n")
for (i in seq_along(percentiles)) {
  thr <- quantile(imp_prefilt, percentiles[i])
  sel <- names(imp_prefilt)[imp_prefilt >= thr]
  rf <- ranger(x = geno[, sel, drop = FALSE], y = y,
               num.trees = ntree, mtry = floor(sqrt(length(sel))),
               importance = "none", oob.error = TRUE)
  subset_results[i, "NLokus"] <- length(sel)
  subset_results[i, "OOB"] <- rf$prediction.error
  write.table(subset_results[i, ], file = "C:/Users/Pc/Desktop/randomf/ranger/OUTPUT/subset_results.tsv",
              sep = "\t", quote = FALSE, row.names = FALSE,
              col.names = FALSE, append = TRUE)
  cat(sprintf("   Percentile %.2f: %d SNP, OOB=%.4f\n", percentiles[i], length(sel), rf$prediction.error))
}

# The best set
best_idx <- which.min(subset_results$OOB)
best_thr <- quantile(imp_prefilt, percentiles[best_idx])
best_snps <- names(imp_prefilt)[imp_prefilt >= best_thr]
write.table(best_snps, file = "C:/Users/Pc/Desktop/randomf/ranger/OUTPUT/best_subset_snps.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# Backward purging
cat("[5/6] Backward purging is starting...\n")
n_total <- length(best_snps)
names_all_iterations <- vector("list", n_total)
names_all_iterations[[n_total]] <- best_snps
write.table(data.frame(Iteration = seq(n_total, 2, by = -1), NLokus = NA, OOB = NA),
            file = "C:/Users/Pc/Desktop/randomf/ranger/OUTPUT/purging_progress.tsv", sep = "\t", row.names = FALSE)

purge_snps <- best_snps
for (k in seq(n_total, 2, by = -1)) {
  rf_iter <- ranger(
    x = geno[, purge_snps, drop = FALSE], y = y,
    num.trees = ntree, mtry = floor(sqrt(length(purge_snps))),
    importance = "impurity", oob.error = TRUE)
  oob_err <- rf_iter$prediction.error
  write.table(data.frame(Iteration = k, NLokus = length(purge_snps), OOB = oob_err),
              file = "C:/Users/Pc/Desktop/randomf/ranger/OUTPUT/purging_progress.tsv",
              sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE, append = TRUE)
  drop_locus <- names(which.min(rf_iter$variable.importance))
  purge_snps <- setdiff(purge_snps, drop_locus)
  names_all_iterations[[length(purge_snps)]] <- purge_snps
}

# save
purge_df <- read.table("C:/Users/Pc/Desktop/randomf/ranger/OUTPUT/purging_progress.tsv", header = TRUE)
best_it <- purge_df$Iteration[which.min(purge_df$OOB)]
final_snps <- names_all_iterations[[best_it]]
write.table(final_snps, file = "C:/Users/Pc/Desktop/randomf/ranger/OUTPUT/best_purged_snps.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
cat(sprintf("[6/6] Backward purging completed. The best %d SNPs have been saved. (C:/Users/Pc/Desktop/randomf/ranger/OUTPUT/best_purged_snps.txt).\n", length(final_snps)))
