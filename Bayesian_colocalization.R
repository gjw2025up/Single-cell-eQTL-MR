###Bayesian colocalization analysis (taking astrocytes as an example)###

# Load required libraries
library(tidyverse)
library(remotes)
library(dplyr)
library(data.table)
library(R.utils)
library(TwoSampleMR)
library(coloc)

# Function to convert effect allele frequency (EAF) to minor allele frequency (MAF)
eaf2maf <- function(eaf = NULL) {
  if (any(is.infinite(eaf))) {
    warning("The 'eaf' vector contains infinite values. These will be converted to NAs.")
    is.na(eaf) <- is.infinite(eaf)
  }
  if (is.null(eaf)) {
    stop("No 'eaf' vector provided.")
  }
  if (is.character(eaf)) {
    stop("'eaf' must be a numeric vector.")
  }
  if (all(is.na(eaf))) {
    stop("All values in 'eaf' are NA.")
  }
  if (is.numeric(eaf)) {
    maf <- eaf
    ind <- which(eaf > 0.5)
    maf[ind] <- 1 - eaf[ind]
    return(maf)
  }
}

# Read raw single-cell eQTL data for astrocytes
ast_eqtl <- fread("celltype-eqtl-sumstats.Ast.tsv.gz", header = T)

# Load and filter MR results to retain significant genes
mr_gene <- read.csv("MR_ast_eqtl_ADHD.csv", header = TRUE)
p_threshold <- 0.05 / length(unique(mr_gene$exposure))  # Bonferroni correction for multiple testing
mr_gene <- subset(mr_gene, pval < p_threshold)
mr_gene <- data.frame(index = unique(mr_gene$exposure))

# Merge eQTL data with MR-selected genes based on gene symbol
ast_eqtl <- merge(ast_eqtl, mr_gene, by.x = "gene_symbol", by.y = "index")

# Load ADHD GWAS summary statistics (outcome data)
outcome_dat <- fread("ADHD.txt", header = T)

# Prepare input data by merging outcome and eQTL datasets on SNP ID
input <- merge(outcome_dat, ast_eqtl, by.x = "SNP", by.y = "snps")

# Harmonize effect alleles between exposure and outcome datasets
input <- input %>% 
  filter((ALT == effect_allele.outcome & other_allele.outcome == REF) |  # Align effect alleles
           (ALT == other_allele.outcome & REF == effect_allele.outcome))
input <- input %>% 
  mutate(beta.outcome = ifelse(ALT == effect_allele.outcome, beta.outcome, -beta.outcome))  # Adjust beta sign for allele consistency

# Calculate variance for beta estimates and convert EAF to MAF
input$varbeta_outcome <- (input$se.outcome * input$se.outcome)
input$varbeta_exposure <- (input$se * input$se)
input$MAF_exposure <- eaf2maf(input$ALT_AF)

# Initialize output data frame for colocalization results
coloc_sum <- data.frame()

# Perform Bayesian colocalization analysis for each gene
index <- unique(input$gene_symbol)
for (i in index) {
  input1 <- subset(input, input$gene_symbol == i)
  result <- coloc.abf(
    dataset1 = list(
      snp = input1$SNP,
      pvalues = input1$pval.outcome,
      type = "cc",
      beta = input1$beta.outcome,
      varbeta = input1$varbeta_outcome,
      N = input1$samplesize.outcome
    ),
    dataset2 = list(
      snp = input1$SNP,
      pvalues = input1$pvalue,
      beta = input1$beta,
      varbeta = input1$varbeta_exposure,
      MAF = input1$MAF_exposure,
      type = "quant",
      N = 424
    )
  )
  # Format results and append to output
  a <- as.data.frame(result$summary)
  a <- data.frame(t(a))
  a$gene <- i
  coloc_sum <- rbind(coloc_sum, data.frame(a))
}

# Save colocalization results to a CSV file
write.csv(coloc_sum, "Colocalization_ast_ADHD.csv", row.names = FALSE)