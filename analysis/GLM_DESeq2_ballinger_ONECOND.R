library( "DESeq2" )
library(glue)


tissue <- "BAT"
temp <- "cold"

x_path <- "../data/male_"
design_path <- "../data/male_design_oneCond"
results_path <- glue("../results/male_")


# full model with intercepts and interaction terms
x <- t(read.delim(glue("{x_path}{tissue}_{temp}_X.txt"), header = FALSE))
int <- ""

null_hypothesis_vector <- c("","_no_cis","_no_trans")

for (i in 1:3) {
    null <- null_hypothesis_vector[i]
    print(null)
    design <- as.matrix( read.delim( glue("{design_path}{null}.txt"), header=FALSE ) )

    # create fake colData
    num_cols <- ncol(x)
    group <- rep(c("A", "B"), length.out = num_cols)
    colData <- data.frame(group = group)

    # fit
    dds <- DESeqDataSetFromMatrix(countData = x, colData = colData, design = design)
    dds = DESeq(dds,full=design,betaPrior=FALSE)


    coefs_to_save <- coef(dds)
    mu <- assays(dds)[["mu"]]
    dispersions <- dispersions(dds)
    logLik <- rowSums(dnbinom(
        x = x,
        mu = mu,
        size = 1 / dispersions,
        log = TRUE
    ))

    write.csv(coefs_to_save, glue("{results_path}{tissue}_{temp}_deseq2{null}_coefficients.csv"), row.names = FALSE)
    write.csv(mu, glue("{results_path}{tissue}_{temp}_deseq2{null}_fitted_vals.csv"), row.names = FALSE)
    write.csv(logLik, glue("{results_path}{tissue}_{temp}_deseq2{null}_logLik.csv"), row.names = FALSE)
    write.csv(dispersions, glue("{results_path}{tissue}_{temp}_deseq2{null}_dispersions.csv"), row.names = FALSE)
}

