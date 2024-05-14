#Please see Slide 8: of /Lab_Presentations/Lab_Meeting_20240208/Lipka_A_Omnigenic_20240208.pptx
# for a visualization of what I want to do with this simulation



#one.replicate.of.simulation.run <- function(){
  
  #Initialize the parameters for simulating the omnigenic genetic architecture the function
  
  source("Simulating_Omnigenic_Genetic_Architecture/setInputParametersforSimplePhentoypes.R")
  this.input.SNPs <- hapmap.file.of.founder.QTLs
  
  these.trait.var.covar <- list()
  these.genetic.value.var.covar <- list()
  these.breeding.value.var.covar <- list()
  
  #Run the function that will conduct the simulations
  source("Four_trait_trick_code/Four_trait_trick_Master_Script_20240509.R")
  
 
  
  #Put the simulated phentoypes from simplePHNEOTYPES back into AlphaSimR
  TrainPop@pheno <- as.matrix(four.traits.omni.core.peri.coreperi$this.simulated.trait[,2])
    
  #Add the population variance-covariance matrix of trait.values, genetic.values, and breeding values to the list
  these.trait.var.covar[[1]] <- var.covar.of.trait.values
  names(these.trait.var.covar)[1] <- "Founder.pop"
  
  these.genetic.value.var.covar[[1]] <- var.covar.of.genetic.values
  names(these.genetic.value.var.covar)[1] <- "Founder.pop"
  
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
    SP.within.function = SP,
    trait.var.covar = these.trait.var.covar,
    genetic.value.var.covar = these.genetic.value.var.covar,
    initial.count.value = length(these.genetic.value.var.covar),
  )
   
  #Join the variance-covariance matrices together
  these.trait.var.covar <- c(these.trait.var.covar, directional.selection.population$trait.var.covar)
  these.genetic.value.var.covar <- c(these.genetic.value.var.covar, 
                                     directional.selection.population$genetic.value.var.covar)
  
  
  #Simulate this.nGenerations generations of disruptive selection
  disruptive.selection.population <- cross.stuff.for.a.whole.bunch.of.generations(
    current.generation = TrainPop,
    type.of.selection = "Disruptive",
    nSelect = this.nSelect,
    nCross = this.nCross,
    nGenerations = this.nGenerations,
    nQtl = this.nQtl,
    SP.within.function = SP,
    trait.var.covar = these.trait.var.covar,
    genetic.value.var.covar = these.genetic.value.var.covar,
    initial.count.value = length(these.genetic.value.var.covar)
  )
  
  #Join the variance-covariance matrices together
  these.trait.var.covar <- c(these.trait.var.covar, disruptive.selection.population$trait.var.covar)
  these.genetic.value.var.covar <- c(these.genetic.value.var.covar, 
                                     disruptive.selection.population$genetic.value.var.covar)
  
  
  #Simulate this.nGenerations generations of stabilizing selection
  stabilizing.selection.population <- cross.stuff.for.a.whole.bunch.of.generations(
    current.generation = TrainPop,
    type.of.selection = "Stabilizing",
    nSelect = this.nSelect,
    nCross = this.nCross,
    nGenerations = this.nGenerations,
    nQtl = this.nQtl,
    SP.within.function = SP,
    trait.var.covar = these.trait.var.covar,
    genetic.value.var.covar = these.genetic.value.var.covar,
    initial.count.value = length(these.genetic.value.var.covar)
  )
  #Join the variance-covariance matrices together
  these.trait.var.covar <- c(these.trait.var.covar,  stabilizing.selection.population$trait.var.covar)
  these.genetic.value.var.covar <- c(these.genetic.value.var.covar, 
                                     stabilizing.selection.population$genetic.value.var.covar)
  
  
  
  
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
  save.image(paste(factor.A,".FactorA..",factor.B, ".FactorB..",
                   factor.C,".FactorC..",factor.D,".FactorD..",
                   this.rep, ".Rep.Rdata", sep = ""))
  
#}#End one.replicate.of.simulation.run
