


#Set your working directory
setwd("/Users/alipka/Library/CloudStorage/Box-Box/Sabbatical_Roslin_Institute/R_workspace/Sabbatical_Project")
home.dir <- getwd()
############################################################
# This code below was written by Chris Gaynor on
# April 27, 2024
#####Read in all of the packages that are necessary
library(AlphaSimR)

# Create founder haplotypes for 100 inbred lines
# Using 2001 loci per chromosome
#   1000 loci for SNPs
#   1 locus for a core gene (i.e. 1 per chromosome)
#   1000 loci for peripheral genes
founderPop = runMacs(100, 10, 2001, inbred=TRUE)

SP = SimParam$new(founderPop)

# Adding SNP chip, these loci will be block from serving as QTL
SP$addSnpChip(1000)

# Adding core genes with a target variance of 0.3 (30% of 1)
SP$addTraitA(1, var=0.3)

# Blocking core genes from being selected as peripheral genes
SP$restrSegSites(excludeQtl = getQtlMap()$id)

# Add 10 traits with epistasis to serve as peripheral gene networks
# Aiming for a total variance of approximately 0.7
for(i in seq_len(10)){
  SP$addTraitAE(1000, relAA=0.3, var=0.7/10)
}

# Add residual error to get a target heritability of 0.5
SP$setVarE(h2 = rep(0.5, 11))

# Creeate initial populations
pop = newPop(founderPop)

# Look at variances

# Additive genetic variance (approximately 1)
popVar(bv(pop) %*% rep(1, 11))

# Total genetic variance
popVar(gv(pop) %*% rep(1, 11))

# Phenotypic variance
popVar(pheno(pop) %*% rep(1, 11))
