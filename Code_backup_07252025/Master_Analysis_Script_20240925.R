#Please see Slide 8: of /Lab_Presentations/Lab_Meeting_20240208/Lipka_A_Omnigenic_20240208.pptx
# for a visualization of what I want to do with this simulation


######I need added code here because I also need QTN and SNP information for the previous generations
directional.subpopulation.QTNs.10.prev.gen <- pullQtlGeno(pop=directional.selection.population.10.pct$previous.generation, 
                                                             trait = 1, asRaw = FALSE, 
                                                             simParam = SP)
directional.subpopulation.QTNs.20.prev.gen <- pullQtlGeno(pop=directional.selection.population.20.pct$previous.generation, 
                                                          trait = 1, asRaw = FALSE, 
                                                          simParam = SP)
disruptive.subpopulation.QTNs.10.prev.gen <- pullQtlGeno(pop=disruptive.selection.population.10.pct$previous.generation, 
                                                         trait = 1, asRaw = FALSE, 
                                                         simParam = SP)
disruptive.subpopulation.QTNs.20.prev.gen <- pullQtlGeno(pop=disruptive.selection.population.20.pct$previous.generation, 
                                                         trait = 1, asRaw = FALSE, 
                                                         simParam = SP) 
stabilizing.subpopulation.QTNs.10.prev.gen <- pullQtlGeno(pop=stabilizing.selection.population.10.pct$previous.generation, 
                                                          trait = 1, asRaw = FALSE, 
                                                          simParam = SP)
stabilizing.subpopulation.QTNs.20.prev.gen <- pullQtlGeno(pop=stabilizing.selection.population.20.pct$previous.generation, 
                                                          trait = 1, asRaw = FALSE, 
                                                          simParam = SP)



directional.subpopulation.10.prev.gen.SNPs  <- pullSnpGeno(directional.selection.population.10.pct$previous.generation, 
                                                           simParam = SP)
directional.subpopulation.20.prev.gen.SNPs <- pullSnpGeno(directional.selection.population.20.pct$previous.generation, 
                                                          simParam = SP)
disruptive.subpopulation.10.prev.gen.SNPs <- pullSnpGeno(disruptive.selection.population.10.pct$previous.generation, 
                                                         simParam = SP)
disruptive.subpopulation.20.prev.gen.SNPs <- pullSnpGeno(disruptive.selection.population.20.pct$previous.generation, 
                                                        simParam = SP)
stabilizing.subpopulation.10.prev.gen.SNPs <- pullSnpGeno(stabilizing.selection.population.10.pct$previous.generation, 
                                                          simParam = SP)
stabilizing.subpopulation.20.prev.gen.SNPs <- pullSnpGeno(stabilizing.selection.population.20.pct$previous.generation, 
                                                          simParam = SP)




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


