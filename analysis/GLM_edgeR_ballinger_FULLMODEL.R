library(edgeR)
library(glue)

x_path <- "../data/ballinger/male_"
design_path <- "../data/ballinger/male_design_model"
results_path <- "../results/ballinger/male_"


# full model with intercepts and interaction terms
tissue <- "LiverBAT"
temp <- "warmcold"
x <- t(read.delim(glue("{x_path}{tissue}_{temp}_X.txt"), header = FALSE))
y <- DGEList(counts=x)  # or the quasi-fitlikelihood F-test with glmQLFit 
y <- normLibSizes(y)
y <- estimateDisp(y)
int <- "int"

null_hypothesis_vector <- c("","_no_cis","_no_trans", 
                            "_no_organ", "_no_temp",
                            "_no_organ-cis","_no_organ-trans",
                            "_no_temp-cis","_no_temp-trans",
                            "_no_organ-temp",
                            "_no_organ-temp-cis","_no_organ-temp-trans")


for (i in 1:12) {
    null <- null_hypothesis_vector[i]
    print(null)
    design <- read.delim( glue("{design_path}_tANDt_{int}_FULL{null}.txt"), header=FALSE ) 
    fit <- glmFit(y,design) # or and glmQLFtest

    coefs_to_save <- coef(fit)
    fitted_vals <- fit$fitted.values
    dev <- fit$deviance
    logLik_dev <- -fit$deviance / 2 
    dispersions <- y$tagwise.dispersion
    logLik <- rowSums(dnbinom(
        x = x,
        mu = fitted_vals,
        size = 1/ dispersions,
        log = TRUE 
    ))


    write.csv(coefs_to_save, glue("{results_path}{tissue}_{temp}_edgeR_{int}_FULL{null}_coefficients.csv"), row.names = FALSE)
    write.csv(fitted_vals, glue("{results_path}{tissue}_{temp}_edgeR_{int}_FULL{null}_fitted_vals.csv"), row.names = FALSE)
    write.csv(logLik, glue("{results_path}{tissue}_{temp}_edgeR_{int}_FULL{null}_logLik.csv"), row.names = FALSE)
    write.csv(logLik_dev, glue("{results_path}{tissue}_{temp}_edgeR_{int}_FULL{null}_logLikdev.csv"), row.names = FALSE)
    write.csv(logLik_dev, glue("{results_path}{tissue}_{temp}_edgeR_{int}_FULL{null}_dev.csv"), row.names = FALSE)
    write.csv(dispersions, glue("{results_path}{tissue}_{temp}_edgeR_{int}_FULL{null}_dispersions.csv"), row.names = FALSE)
}