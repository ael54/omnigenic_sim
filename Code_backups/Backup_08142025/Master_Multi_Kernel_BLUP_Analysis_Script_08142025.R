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
###Genomic prediction across subpopulations phase

#For the QTLs
source("QG_Analysis_Code/Run_Multi_Kernel_BLUP_GS_on_QTNs.R")

#For the SNPs
source("QG_Analysis_Code/Run_Multi_Kernel_BLUP_GS_on_SNPs.R")


