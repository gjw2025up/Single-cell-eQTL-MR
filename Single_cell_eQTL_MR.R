####Single-cell eQTL MR analysis (taking astrocytes as an example)####

# Load necessary libraries
library(tidyverse)
library(data.table)
library(R.utils)
library(devtools)
library(ieugwasr)
library(TwoSampleMR)
library(doParallel)

# Set the OpenGWAS JWT token. Replace "Your token" with the actual token.
# https://api.opengwas.io/api/
Sys.setenv(OPENGWAS_JWT = "Your token")
get_opengwas_jwt()
user()

# Extract eQTL data and filter by FDR
ast_eqtl <- fread("celltype-eqtl-sumstats.Ast.tsv.gz", header = T)
ast_eqtl <- subset(ast_eqtl, ast_eqtl$significant_by_2step_FDR == "Yes")

# Select relevant columns and rename them
ast_eqtl <- ast_eqtl[, c("snps", "ALT", "REF", "ALT_AF", "beta", "se", "pvalue", "celltype", "gene_symbol", "gene_id")]
colnames(ast_eqtl) <- c("SNP", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure", "beta.exposure", "se.exposure", "pval.exposure", "celltype", "exposure", "gene_id")
ast_eqtl$id.exposure <- ast_eqtl$gene_id
ast_eqtl$samplesize.exposure <- 424

# Perform LD clumping
exposure_dat <- clump_data(ast_eqtl, clump_kb = 10000,
                           clump_r2 = 0.1,
                           clump_p1 = 1,
                           clump_p2 = 1,
                           pop = "EUR")

# write.csv(dat,file = "ast_eqtl_fdr0.05_clump_10000kb_0.1r2.csv")

# Preprocess outcome data (using ADHD GWAS iPSYCH_deCODE_PGC as an example)
outcome_dat <- fread("ADHD2022_iPSYCH_deCODE_PGC.meta.gz", header = T)
outcome_dat$beta = log(outcome_dat$OR) # Calculate beta from OR
outcome_dat$N <- outcome_dat$Nca + outcome_dat$Nco # Calculate the total sample size
# Select relevant columns and rename them
outcome_dat <- outcome_dat[, c("SNP", "CHR", "BP", "A1", "A2", "FRQ_A_38691", "N", "beta", "SE", "P")]
colnames(outcome_dat) <- c("SNP", "chr", "pos", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "samplesize.outcome", "beta.outcome", "se.outcome", "pval.outcome")
outcome_dat$outcome <- "ADHD"
outcome_dat$id.outcome <- "ADHD"

# write.table(dat,"ADHD.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

# Harmonize exposure and outcome data
dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)

# Calculate R2 and F statistics
dat$R2 <- (2 * dat$beta.exposure * dat$beta.exposure * dat$eaf.exposure * (1 - dat$eaf.exposure) / (2 * dat$beta.exposure * dat$beta.exposure * dat$eaf.exposure * (1 - dat$eaf.exposure) + 2 * dat$se.exposure * dat$se.exposure * dat$samplesize.exposure * dat$eaf.exposure * (1 - dat$eaf.exposure)))
dat$F <- dat$R2 * (dat$samplesize.exposure - 2) / (1 - dat$R2)

# Perform Steiger filtering
steiger_dat <- steiger_filtering(dat)
dat <- subset(steiger_dat, steiger_dir == TRUE)

# Exclude SNPs significantly associated with the outcome
dat <- subset(dat,pval.outcome >= 5e-8)

# Split the data by exposure
list_dat <- split(dat, dat$exposure)

# Define a function for MR analysis
mr_Speed <- function(dat) {
  res = TwoSampleMR::mr(
    dat,
    method_list = c("mr_egger_regression",
                    "mr_weighted_median",
                    "mr_ivw",
                    "mr_ivw_mre",
                    "mr_ivw_fe",
                    "mr_wald_ratio")
  )
  return(res)
}


detectCores(logical = FALSE) # Check the number of available CPU cores
cl <- makeCluster(5) # Create a cluster with 5 cores
registerDoParallel(cl) # Register the cluster for parallel processing
clusterExport(cl = cl, varlist = ls()) # Export all variables to the cluster
start1 = Sys.time() # Record the start time

# Perform parallel MR analysis
res = parLapply(cl = cl,
                X = list_dat,
                fun = mr_Speed)
end1 = Sys.time(); print(end1 - start1) # Record the end time and print the elapsed time

res_df <- do.call(rbind, res) # Combine the results into a single data frame

# Write the results to a CSV file
write.csv(res_df, file = "MR_ast_eqtl_ADHD.csv")