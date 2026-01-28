library(devtools)
library(edgeR)
devtools::load_all()



metadata_path_full <- system.file("extdata", "ballinger_metadata.csv", package = "XgeneR")
metadata_full <- read.csv(metadata_path_full, row.names = 1) 
rownames(metadata_full) <- metadata_full$Sample

count_path_full <- system.file("extdata", "cold_ballinger_counts.csv", package = "XgeneR")
count_path_cold <- system.file("extdata", "cold_ballinger_counts.csv", package = "XgeneR")
count_path_warm <- system.file("extdata", "warm_ballinger_counts.csv", package = "XgeneR")


counts_full <- read.csv(count_path_full, row.names = 1)
counts_full <- as.matrix(counts_full)
counts_cold <- read.csv(count_path_cold, row.names = 1)
counts_cold <- as.matrix(counts_cold)
counts_warm <- read.csv(count_path_warm, row.names = 1)
counts_warm <- as.matrix(counts_warm)

# Cold, free model
metadata_ <- metadata_full[metadata_full$Temperature == "Cold", , drop = FALSE]

fit_obj <- new("fitObject",counts = counts_cold, metadata = metadata_, covariate_cols = c("Tissue","Individual"))
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
write.csv(fit_obj@weights, file = "../results/cold_weights_free.csv",row.names = TRUE)
write.csv(out_df,file = "../results/cold_sig_free.csv",row.names = TRUE)
write.csv(fit_obj@design_matrix_full, file = "../results/cold_design_matrix_free.csv",row.names = TRUE)

# Cold, log-additive model
fit_obj <- new("fitObject",counts = counts_cold, metadata = metadata_, covariate_cols = c("Tissue"))
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
write.csv(fit_obj@weights, file = "../results/cold_weights_log_additive.csv",row.names = TRUE)
write.csv(out_df,file = "../results/cold_sig_log_additive.csv",row.names = TRUE)
write.csv(fit_obj@design_matrix_full, file = "../results/cold_design_matrix_log_additive.csv",row.names = TRUE)


# # Cold, dominant model
# fit_obj <- new("fitObject",counts = counts_cold, metadata = metadata_, covariate_cols = c("Tissue"), trans_model = "dominant")
# fit_obj <- fit_edgeR(fit_obj)

# weight_names <- setdiff(colnames(fit_obj@raw_pvals), "Genes")

# out_df <- data.frame(
#   setNames(
#     fit_obj@raw_pvals[, weight_names, drop = FALSE],
#     paste0("raw_pval_", weight_names)
#   ),
#   setNames(
#     fit_obj@BH_FDRs[, weight_names, drop = FALSE],
#     paste0("BH_FDR_", weight_names)
#   ),
#   row.names = rownames(fit_obj@raw_pvals)
# )

# # write results
# write.csv(fit_obj@weights, file = "/home/maria/igvf/PBMC/results/ballinger/cold_weights_dominant.csv",row.names = TRUE)
# write.csv(out_df,file = "/home/maria/igvf/PBMC/results/ballinger/cold_sig_dominant.csv",row.names = TRUE)
# write.csv(fit_obj@design_matrix_full, file = "/home/maria/igvf/PBMC/results/ballinger/cold_design_matrix_dominant.csv",row.names = TRUE)



# Cold, BAT, log-additive 
metadata_ <- metadata_full[metadata_full$Temperature == "Cold" & metadata_full$Tissue == "BAT", , drop = FALSE]
counts_ <- counts_cold[, rownames(metadata_), drop = FALSE]


fit_obj <- new("fitObject",counts = counts_, metadata = metadata_)
fit_obj <- fit_edgeR(fit_obj)

design <- fit_obj@design_matrix_full
y <- DGEList(counts=counts_)  # or the quasi-fitlikelihood F-test with glmQLFit 
y <- normLibSizes(y)
y <- estimateDisp(y)

fit <- glmFit(y,design) # or and glmQLFtest
fitted_vals <- fit$fitted.values
logLik_dev <- -fit$deviance / 2 
dispersions <- y$tagwise.dispersion
logLik <- rowSums(dnbinom(
    x = counts_,
    mu = fitted_vals,
    size = 1/ dispersions,
    log = TRUE
))

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

out_df$logLik <- logLik
out_df$dispersions <- dispersions
 


write.csv(fitted_vals, "../results/BATcold_fitted_vals_log_additive.csv", row.names = FALSE)
write.csv(fit_obj@weights, file = "../results/BATcold_weights_log_additive.csv",row.names = TRUE)
write.csv(out_df,file = "../results/BATcold_sig_log_additive.csv",row.names = TRUE)
write.csv(fit_obj@design_matrix_full, file = "../results/BATcold_design_matrix_log_additive.csv",row.names = TRUE)


# Cold, BAT, Free 
metadata_ <- metadata_full[metadata_full$Temperature == "Cold" & metadata_full$Tissue == "BAT", , drop = FALSE]
counts_ <- counts_cold[, rownames(metadata_), drop = FALSE]


fit_obj <- new("fitObject",counts = counts_, metadata = metadata_, covariate_cols = c("Individual"))
fit_obj <- fit_edgeR(fit_obj)

design <- fit_obj@design_matrix_full
y <- DGEList(counts=counts_)  # or the quasi-fitlikelihood F-test with glmQLFit 
y <- normLibSizes(y)
y <- estimateDisp(y)

fit <- glmFit(y,design) # or and glmQLFtest
fitted_vals <- fit$fitted.values
logLik_dev <- -fit$deviance / 2 
dispersions <- y$tagwise.dispersion
logLik <- rowSums(dnbinom(
    x = counts_,
    mu = fitted_vals,
    size = 1/ dispersions,
    log = TRUE
))

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

out_df$logLik <- logLik
out_df$dispersions <- dispersions
 


write.csv(fitted_vals, "../results/BATcold_fitted_vals_free.csv", row.names = FALSE)
write.csv(fit_obj@weights, file = "../resultsBATcold_weights_free.csv",row.names = TRUE)
write.csv(out_df,file = "../results/BATcold_sig_free.csv",row.names = TRUE)
write.csv(fit_obj@design_matrix_full, file = "../results/BATcold_design_matrix_free.csv",row.names = TRUE)

# Cold, BAT, dominant
metadata_ <- metadata_full[metadata_full$Temperature == "Cold" & metadata_full$Tissue == "BAT", , drop = FALSE]
counts_ <- counts_cold[, rownames(metadata_), drop = FALSE]


fit_obj <- new("fitObject",counts = counts_, metadata = metadata_, trans_model = "dominant")
fit_obj <- fit_edgeR(fit_obj)

design <- fit_obj@design_matrix_full
y <- DGEList(counts=counts_)  # or the quasi-fitlikelihood F-test with glmQLFit 
y <- normLibSizes(y)
y <- estimateDisp(y)

fit <- glmFit(y,design) # or and glmQLFtest
fitted_vals <- fit$fitted.values
logLik_dev <- -fit$deviance / 2 
dispersions <- y$tagwise.dispersion
logLik <- rowSums(dnbinom(
    x = counts_,
    mu = fitted_vals,
    size = 1/ dispersions,
    log = TRUE
))

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

out_df$logLik <- logLik
out_df$dispersions <- dispersions
 


write.csv(fitted_vals, "../results/BATcold_fitted_vals_dominant.csv", row.names = FALSE)
write.csv(fit_obj@weights, file = "../results/ballinger/BATcold_weights_dominant.csv",row.names = TRUE)
write.csv(out_df,file = "../results/BATcold_sig_dominant.csv",row.names = TRUE)
write.csv(fit_obj@design_matrix_full, file = "../results/BATcold_design_matrix_dominant.csv",row.names = TRUE)


# Cold, Liver, log-additive 
metadata_ <- metadata_full[metadata_full$Temperature == "Cold" & metadata_full$Tissue == "Liver", , drop = FALSE]
counts_ <- counts_cold[, rownames(metadata_), drop = FALSE]


fit_obj <- new("fitObject",counts = counts_, metadata = metadata_)
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

write.csv(fit_obj@weights, file = "../results/Livercold_weights_log_additive.csv",row.names = TRUE)
write.csv(out_df,file = "../results/Livercold_sig_log_additive.csv",row.names = TRUE)
write.csv(fit_obj@design_matrix_full, file = "../results/Livercold_design_matrix_log_additive.csv",row.names = TRUE)

# Warm, BAT, log-additive 
metadata_ <- metadata_full[metadata_full$Temperature == "Warm" & metadata_full$Tissue == "BAT", , drop = FALSE]
counts_ <- counts_warm[, rownames(metadata_), drop = FALSE]


fit_obj <- new("fitObject",counts = counts_, metadata = metadata_)
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

write.csv(fit_obj@weights, file = "../results/BATwarm_weights_log_additive.csv",row.names = TRUE)
write.csv(out_df,file = "../resultsBATwarm_sig_log_additive.csv",row.names = TRUE)
write.csv(fit_obj@design_matrix_full, file = "../results/BATwarm_design_matrix_log_additive.csv",row.names = TRUE)

# Warm, Liver, log-additive 
metadata_ <- metadata_full[metadata_full$Temperature == "Warm" & metadata_full$Tissue == "Liver", , drop = FALSE]
counts_ <- counts_warm[, rownames(metadata_), drop = FALSE]


fit_obj <- new("fitObject",counts = counts_, metadata = metadata_)
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

write.csv(fit_obj@weights, file = "../results/Liverwarm_weights_log_additive.csv",row.names = TRUE)
write.csv(out_df,file = "../results/Liverwarm_sig_log_additive.csv",row.names = TRUE)
write.csv(fit_obj@design_matrix_full, file = "../results/Liverwarm_design_matrix_log_additive.csv",row.names = TRUE)
