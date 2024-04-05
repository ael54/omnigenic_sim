#Please see Slide 8: of /Lab_Presentations/Lab_Meeting_20240208/Lipka_A_Omnigenic_20240208.pptx
# for a visualization of what I want to do with this simualtion

#Clear your workspace
rm(list = ls())

#Set the working directory
setwd("/Users/alipka/Library/CloudStorage/Box-Box/Sabbatical_Roslin_Institute/R_workspace/Sabbatical_Project")
home.dir <- getwd()



###Read in the R workspace of simulated data
load("2.FactorA..0.05.FactorB..0.05.FactorC..0.05.FactorD..3.Rep.Rdata")


#####Read in all of the packages that are necessary
#Read in prerequiste libaries for GAPIT
library('MASS')
library(multtest)
library(gplots)
library(rrBLUP)



#Read in GAPIT
setwd("Scripts_Necessary_for_GAPIT")
source("GAPIT_EMMA source code.txt")
source("GAPIT_Code_from_Internet_20120411_Allelic_Effect.r")
setwd(home.dir)



################################################################################
################################################################################
###GWAS phase
#For the QTNs

#Calculate the spearman correlation coefficients between additive 
# effect estimates core and peripheral QTNs
source("QG_Analysis_Code/Run_GWAS_on_QTNs.R")

#Object for spearman correlation coefficients between core and peripheral QTNs
spearman.correlations.between.core.QTNs
spearman.correlations.between.peripheral.QTNs


#Calculate the spearman correlation coefficients between additive 
# effect estimates of "Core SNPs" that are within 0.05 cM of core QTNs
# and "Peripheral SNPs that are not"
source("QG_Analysis_Code/Run_GWAS_on_SNPs.R")

#Object for spearman correlation coefficients between core and peripheral QTNs
spearman.correlations.between.core.SNPs
spearman.correlations.between.peripheral.SNPs

 
################################################################################
################################################################################
###Estimating variance components phase

#For the QTLs




#For the SNPs





################################################################################
################################################################################
###Genomic prediction across subpopulations phase

#For the QTLs
source("QG_Analysis_Code/Run_GS_on_QTNs.R")



#For the SNPs
source("QG_Analysis_Code/Run_GS_on_SNPs.R")


