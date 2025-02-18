library(edgeR)
library(glue)

x_path <- "../data/ballinger/male_"
design_path <- "../data/ballinger/male_design_oneCond.txt"
results_path <- "../results/ballinger/male_"

# design matrix
design <- read.delim( design_path, header=FALSE )

# Male BAT cold
x <- t(read.delim(glue("{x_path}BAT_cold_X.txt"), header = FALSE))
y <- DGEList(counts=x)  # or the quasi-fitlikelihood F-test with glmQLFit 
y <- normLibSizes(y)
y <- estimateDisp(y)


fit <- glmFit(y,design) # or and glmQLFtest
lrt_no_cis <- glmLRT(fit,coef=2)
lrt_no_trans <- glmLRT(fit,contrast=c(0,0,-1,1))
coefs_to_save <- coef(fit)
fitted_vals <- fit$fitted.values
logLik_dev <- -fit$deviance / 2 
dispersions <- y$tagwise.dispersion
logLik <- rowSums(dnbinom(
    x = x,
    mu = fitted_vals,
    size = 1/ dispersions,
    log = TRUE
))

pvals_no_cis <- lrt_no_cis$table$PValue
pvals_no_trans <- lrt_no_trans$table$PValue

write.csv(logLik, glue("{results_path}BAT_cold_edgeR_logLik.csv"), row.names = FALSE)
write.csv(logLik_dev, glue("{results_path}BAT_cold_edgeR_logLikdev.csv"), row.names = FALSE)
write.csv(pvals_no_cis, glue("{results_path}BAT_cold_edgeR_pvals_no_cis.csv"), row.names = FALSE)
write.csv(pvals_no_trans, glue("{results_path}BAT_cold_edgeR_pvals_no_trans.csv"), row.names = FALSE)
write.csv(coefs_to_save, glue("{results_path}BAT_cold_edgeR_coefficients.csv"), row.names = FALSE)
write.csv(fitted_vals, glue("{results_path}BAT_cold_edgeR_fitted_vals.csv"), row.names = FALSE)
write.csv(dispersions, glue("{results_path}BAT_cold_edgeR_dispersions.csv"), row.names = FALSE)

# Male BAT warm
x <- t(read.delim(glue("{x_path}BAT_warm_X.txt"), header = FALSE))
y <- DGEList(counts=x)  # or the quasi-fitlikelihood F-test with glmQLFit 
y <- normLibSizes(y)
y <- estimateDisp(y)


fit <- glmFit(y,design) # or and glmQLFtest
lrt_no_cis <- glmLRT(fit,coef=2)
lrt_no_trans <- glmLRT(fit,contrast=c(0,0,-1,1))
coefs_to_save <- coef(fit)
fitted_vals <- fit$fitted.values
logLik_dev <- -fit$deviance / 2 
dispersions <- y$tagwise.dispersion
logLik <- rowSums(dnbinom(
    x = x,
    mu = fitted_vals,
    size = 1/ dispersions,
    log = TRUE
))

pvals_no_cis <- lrt_no_cis$table$PValue
pvals_no_trans <- lrt_no_trans$table$PValue

write.csv(logLik, glue("{results_path}BAT_warm_edgeR_logLik.csv"), row.names = FALSE)
write.csv(logLik_dev, glue("{results_path}BAT_warm_edgeR_logLikdev.csv"), row.names = FALSE)
write.csv(pvals_no_cis, glue("{results_path}BAT_warm_edgeR_pvals_no_cis.csv"), row.names = FALSE)
write.csv(pvals_no_trans, glue("{results_path}BAT_warm_edgeR_pvals_no_trans.csv"), row.names = FALSE)
write.csv(coefs_to_save, glue("{results_path}BAT_warm_edgeR_coefficients.csv"), row.names = FALSE)
write.csv(fitted_vals, glue("{results_path}BAT_warm_edgeR_fitted_vals.csv"), row.names = FALSE)

# Male Liver cold
x <- t(read.delim(glue("{x_path}Liver_cold_X.txt"), header = FALSE))
y <- DGEList(counts=x)  # or the quasi-fitlikelihood F-test with glmQLFit 
y <- normLibSizes(y)
y <- estimateDisp(y)


fit <- glmFit(y,design) # or and glmQLFtest
lrt_no_cis <- glmLRT(fit,coef=2)
lrt_no_trans <- glmLRT(fit,contrast=c(0,0,-1,1))
coefs_to_save <- coef(fit)
fitted_vals <- fit$fitted.values
logLik_dev <- -fit$deviance / 2 
dispersions <- y$tagwise.dispersion
logLik <- rowSums(dnbinom(
    x = x,
    mu = fitted_vals,
    size = 1/ dispersions,
    log = TRUE
))

pvals_no_cis <- lrt_no_cis$table$PValue
pvals_no_trans <- lrt_no_trans$table$PValue

write.csv(logLik, glue("{results_path}Liver_cold_edgeR_logLik.csv"), row.names = FALSE)
write.csv(logLik_dev, glue("{results_path}Liver_cold_edgeR_logLikdev.csv"), row.names = FALSE)
write.csv(pvals_no_cis, glue("{results_path}Liver_cold_edgeR_pvals_no_cis.csv"), row.names = FALSE)
write.csv(pvals_no_trans, glue("{results_path}Liver_cold_edgeR_pvals_no_trans.csv"), row.names = FALSE)
write.csv(coefs_to_save, glue("{results_path}Liver_cold_edgeR_coefficients.csv"), row.names = FALSE)
write.csv(fitted_vals, glue("{results_path}Liver_cold_edgeR_fitted_vals.csv"), row.names = FALSE)


# Male Liver warm
x <- t(read.delim(glue("{x_path}Liver_warm_X.txt"), header = FALSE))
y <- DGEList(counts=x)  # or the quasi-fitlikelihood F-test with glmQLFit 
y <- normLibSizes(y)
y <- estimateDisp(y)


fit <- glmFit(y,design) # or and glmQLFtest
lrt_no_cis <- glmLRT(fit,coef=2)
lrt_no_trans <- glmLRT(fit,contrast=c(0,0,-1,1))
coefs_to_save <- coef(fit)
fitted_vals <- fit$fitted.values
logLik_dev <- -fit$deviance / 2 
dev <- fit$deviance
dispersions <- y$tagwise.dispersion
logLik <- rowSums(dnbinom(
    x = x,
    mu = fitted_vals,
    size = 1/ dispersions,
    log = TRUE
))

pvals_no_cis <- lrt_no_cis$table$PValue
pvals_no_trans <- lrt_no_trans$table$PValue

write.csv(logLik, glue("{results_path}Liver_warm_edgeR_dev.csv"), row.names = FALSE)
write.csv(logLik, glue("{results_path}Liver_warm_edgeR_logLik.csv"), row.names = FALSE)
write.csv(logLik_dev, glue("{results_path}Liver_warm_edgeR_logLikdev.csv"), row.names = FALSE)
write.csv(pvals_no_cis, glue("{results_path}Liver_warm_edgeR_pvals_no_cis.csv"), row.names = FALSE)
write.csv(pvals_no_trans, glue("{results_path}Liver_warm_edgeR_pvals_no_trans.csv"), row.names = FALSE)
write.csv(coefs_to_save, glue("{results_path}Liver_warm_edgeR_coefficients.csv"), row.names = FALSE)
write.csv(fitted_vals, glue("{results_path}Liver_warm_edgeR_fitted_vals.csv"), row.names = FALSE)

