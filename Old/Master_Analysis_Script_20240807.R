#Please see Slide 8: of /Lab_Presentations/Lab_Meeting_20240208/Lipka_A_Omnigenic_20240208.pptx
# for a visualization of what I want to do with this simualtion




################################################################################
################################################################################
###GWAS phase
#For the QTNs

#Calculate the spearman correlation coefficients between additive 
# effect estimates core and peripheral QTNs
source("QG_Analysis_Code/Run_GWAS_on_QTNs_20240815.R")

#Object for spearman correlation coefficients between core and peripheral QTNs
spearman.correlations.between.core.QTNs
median.spearman.correlations.between.core.QTN
median.core.spearman.correlation.within.types.of.selection.QTN
median.core.spearman.correlation.within.prev.gen.QTN
median.core.spearman.correlation.different.types.of.selection.QTN

spearman.correlations.between.peripheral.QTNs
median.spearman.correlations.between.peripheral.QTN
median.peripheral.spearman.correlation.within.types.of.selection.QTN
median.peripheral.spearman.correlation.within.prev.gen.QTN 
median.peripheral.spearman.correlation.different.types.of.selection.QTN

#Calculate the spearman correlation coefficients between additive 
# effect estimates of "Core SNPs" that are within 0.05 cM of core QTNs
# and "Peripheral SNPs that are not"
source("QG_Analysis_Code/Run_GWAS_on_SNPs_20240815.R")

#Object for spearman correlation coefficients between core and peripheral QTNs
spearman.correlations.between.core.SNPs
median.spearman.correlations.between.core.SNPs
median.core.spearman.correlation.within.types.of.selection.SNPs
median.core.spearman.correlation.within.prev.gen.SNPs
median.core.spearman.correlation.different.types.of.selection.SNPs

spearman.correlations.between.peripheral.SNPs
median.spearman.correlations.between.peripheral.SNPs 
median.peripheral.spearman.correlation.within.types.of.selection.SNPs
median.peripheral.spearman.correlation.within.prev.gen.SNPs
median.peripheral.spearman.correlation.different.types.of.selection.SNPs
 
################################################################################
################################################################################
###Estimating variance components phase

#For the QTLs




#For the SNPs





################################################################################
################################################################################
###Genomic prediction across subpopulations phase

#For the QTLs
source("QG_Analysis_Code/Run_GS_on_QTNs_20240814.R")



#For the SNPs
source("QG_Analysis_Code/Run_GS_on_SNPs_20240814.R")


