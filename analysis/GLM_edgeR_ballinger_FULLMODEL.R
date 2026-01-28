library(devtools)
library(edgeR)
devtools::load_all()


count_path <- system.file("extdata", "ballinger_counts.csv", package = "XgeneR")
metadata_path <- system.file("extdata", "ballinger_metadata.csv", package = "XgeneR")


counts <- read.csv(count_path, row.names = 1)
counts <- as.matrix(counts)
metadata <- read.csv(metadata_path, row.names = 1) 
rownames(metadata) <- metadata$Sample


fit_obj <- new("fitObject",counts = counts, metadata = metadata, fields_to_test = c("Tissue","Temperature"), higher_order_interactions = list("Reg*BAT*Cold","BAT*Cold"))
fit_obj <- fit_edgeR(fit_obj)

weight_names <- setdiff(colnames(fit_obj@raw_pvals), "Genes")

out_df <- data.frame(
  setNames(
    fit_obj@raw_pvals[, weight_names, drop = FALSE],
    paste0("raw_pval_", weight_names)
  ),
  setNames(
    fit_obj@BH_FDRs[, weight_names, drop = FALSE],
    paste0("BH_FDR_", weight_names)
  ),
  row.names = rownames(fit_obj@raw_pvals)
)

# write results
write.csv(fit_obj@weights, file = "../results/weights_log_additive.csv",row.names = TRUE)
write.csv(out_df,file = "../results/sig_log_additive.csv",row.names = TRUE)
write.csv(fit_obj@design_matrix_full, file = "../results/design_matrix_log_additive.csv",row.names = TRUE)