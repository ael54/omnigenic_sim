#Please see Slide 8: of /Lab_Presentations/Lab_Meeting_20240208/Lipka_A_Omnigenic_20240208.pptx
# for a visualization of what I want to do with this simualtion

#Clear your workspace

one.replicate.of.simulation.run <- function(){
  
  #Initialize the parameters for simulating the omnigenic genetic architecture the function
  
  source("Simulating_Omnigenic_Genetic_Architecture/setInputParametersforSimplePhentoypes.R")
  this.input.SNPs <- hapmap.file.of.founder.QTLs
  
  
  #Run the function that will conduct the simulations
  this.simulated.trait <- simulate.omnigenic.architecture(input.SNPs = this.input.SNPs,
                            number.of.trait.reps = this.number.of.trait.reps,
                            broad.sense.H2 = this.broad.sense.H2,
                            number.of.core.genes = this.number.of.core.genes,
                            sigma.2.C = this.sigma.2.C,
                            sigma.2.N = this.sigma.2.N,
                            sigma.2.CC = this.sigma.2.CC,
                            sigma.2.NN = this.sigma.2.NN,
                            sigma.2.CN = this.sigma.2.CN, 
                            seed.number.core.vs.perhiperal = this.seed.number.core.vs.perhiperal,
                            seed.number.core.add = this.seed.number.core.add,
                            seed.number.peripheral.add = this.seed.number.peripheral.add,
                            seed.number.core.core.epi = this.seed.number.core.core.epi,
                            seed.number.peri.peri.epi = this.seed.number.peri.peri.epi,
                            seed.number.core.peri.epi = this.seed.number.core.peri.epi,
                            snps.are.in.columns = this.snps.are.in.columns,
                            seed.within.simplePHENOTYPES = this.rep,
                            output.directory.name = this.output.directory.name)
  
  #Put the simulated phentoypes from simplePHNEOTYPES back into AlphaSimR
  TrainPop@pheno <- as.matrix(this.simulated.trait$this.simulated.trait[,2])
    
  ################################################################################
  ################################################################################
  #For loop through various breeding programs (directional seletion,
  #            disruptive selection, stabilizing selection)... 
  
    
  #Simulate this.nGenerations generations of directional selection
  
  
  directional.selection.population <- cross.stuff.for.a.whole.bunch.of.generations(
    current.generation = TrainPop,
    type.of.selection = "Direct",
    nSelect = this.nSelect,
    nCross = this.nCross,
    nGenerations = this.nGenerations,
    nQtl = this.nQtl,
    SP.within.function = SP
  )
   
  
  #Simulate this.nGenerations generations of disruptive selection
  disruptive.selection.population <- cross.stuff.for.a.whole.bunch.of.generations(
    current.generation = TrainPop,
    type.of.selection = "Disruptive",
    nSelect = this.nSelect,
    nCross = this.nCross,
    nGenerations = this.nGenerations,
    nQtl = this.nQtl,
    SP.within.function = SP
  )
  
  
  #Simulate this.nGenerations generations of stabilizing selection
  stabilizing.selection.population <- cross.stuff.for.a.whole.bunch.of.generations(
    current.generation = TrainPop,
    type.of.selection = "Stabilizing",
    nSelect = this.nSelect,
    nCross = this.nCross,
    nGenerations = this.nGenerations,
    nQtl = this.nQtl,
    SP.within.function = SP
  )
  
  
  
  #Merge the SNP data together from the three different subpopulations
  directional.subpopulation.SNPs <- pullSnpGeno(directional.selection.population, 
                                           simParam = SP)
  
  disruptive.subpopulation.SNPs <- pullSnpGeno(disruptive.selection.population, 
                                           simParam = SP)
  
  stabilizing.subpopulation.SNPs <- pullSnpGeno(stabilizing.selection.population, 
                                               simParam = SP)
  
  combined.subpopulation.SNPs <- rbind(directional.subpopulation.SNPs,
                    disruptive.subpopulation.SNPs,stabilizing.subpopulation.SNPs)
  
  
  #Merge the QTN data together from the three different subpopulations
  directional.subpopulation.QTNs <- pullQtlGeno(pop=directional.selection.population, 
                                                trait = 1, asRaw = FALSE, 
                                                simParam = SP)
  
  
  
  disruptive.subpopulation.QTNs <- pullQtlGeno(pop=disruptive.selection.population, 
                                               trait = 1, asRaw = FALSE, 
                                               simParam = SP)
  
  
  stabilizing.subpopulation.QTNs <- pullQtlGeno(pop=stabilizing.selection.population, 
                                                trait = 1, asRaw = FALSE, 
                                                simParam = SP)
  
  combined.subpopulation.QTNs <- rbind(directional.subpopulation.QTNs,
                    disruptive.subpopulation.QTNs,stabilizing.subpopulation.QTNs)
  
  #Merge the trait data together from the three different subpopulations
  #####The code below is not working...fix it
  directional.subpopulation.trait <- directional.selection.population@pheno
  row.names(directional.subpopulation.trait) = row.names(directional.subpopulation.SNPs)
  
  disruptive.subpopulation.trait <- disruptive.selection.population@pheno
  row.names(disruptive.subpopulation.trait) = row.names(disruptive.subpopulation.SNPs)
  
  
  stabilizing.subpopulation.trait <- stabilizing.selection.population@pheno
  row.names(stabilizing.subpopulation.trait) = row.names(stabilizing.subpopulation.SNPs)
  
  
  combined.subpopulation.trait <- rbind(directional.subpopulation.trait,
                    disruptive.subpopulation.trait,stabilizing.subpopulation.trait)
  
  
  ###Save the R workspace
  save.image("Test_Run_small_data_set_Quick_Haplo_instead_of_runMacs2_20240314.Rdata")

}#End one.replicate.of.simulation.run
