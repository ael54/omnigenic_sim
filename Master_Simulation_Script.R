#Please see Slide 8: of /Lab_Presentations/Lab_Meeting_20240208/Lipka_A_Omnigenic_20240208.pptx
# for a visualization of what I want to do with this simulation



#one.replicate.of.simulation.run <- function(){
  
  #Initialize the parameters for simulating the omnigenic genetic architecture the function
  
  source("Simulating_Omnigenic_Genetic_Architecture/setInputParametersforSimplePhentoypes.R")
  this.input.SNPs <- hapmap.file.of.founder.QTLs
  
  these.trait.var.covar <- list()
  these.genetic.value.var.covar <- list()
  these.breeding.value.var.covar <- list()
  
  these.trait.values <- list()
  these.genetic.values <- list()
  these.breeding.values <- list()
  
  #Run the function that will conduct the simulations
 time.1 <- Sys.time()
  source("Four_trait_trick_code/Four_trait_trick_Master_Script.R")
 time.2 <- Sys.time()
 time.2-time.1
 
  
  #Put the simulated phentoypes from simplePHENOTYPES back into AlphaSimR
  TrainPop@pheno <- as.matrix(four.traits.omni.core.peri.coreperi$this.simulated.trait[,2])
    
  #Add the population variance-covariance matrix of trait.values, genetic.values, and breeding values to the list
  these.trait.var.covar[[1]] <- var.covar.of.trait.values
  names(these.trait.var.covar)[1] <- paste("Founder.pop.Rep.",this.rep, sep = "")
  
  these.genetic.value.var.covar[[1]] <- var.covar.of.genetic.values
  names(these.genetic.value.var.covar)[1] <-  paste("Founder.pop.Rep.",this.rep, sep = "")
  
  these.breeding.value.var.covar[[1]] <- var.covar.of.breeding.values
  names(these.breeding.value.var.covar)[1] <-  paste("Founder.pop.Rep.",this.rep, sep = "")
  
  #Add the observed trait.values, genetic.values, and breeding values to the list 
  these.trait.values[[1]] <- the.trait.values 
  names(these.trait.values)[1] <- paste("Founder.pop.Rep.",this.rep, sep = "")
  
  these.genetic.values[[1]] <- the.genetic.values
  names(these.genetic.values )[1] <-  paste("Founder.pop.Rep.",this.rep, sep = "")
  
  these.breeding.values[[1]] <- the.breeding.values
  names(these.breeding.values)[1] <-  paste("Founder.pop.Rep.",this.rep, sep = "")
  ################################################################################
  ################################################################################
  #For loop through various breeding programs (directional seletion,
  #            disruptive selection, stabilizing selection)... 
  
    
  #Simulate this.nGenerations generations of directional selection, selection on 10% of the individuals
    directional.selection.population.10.pct <- cross.stuff.for.a.whole.bunch.of.generations(
    current.generation = TrainPop,
    type.of.selection = "Direct",
    nSelect = this.nSelect,
    nCross = this.nCross,
    nGenerations = this.nGenerations,
    nQtl = this.nQtl,
    SP.within.function = SP,
    trait.var.covar = these.trait.var.covar,
    genetic.value.var.covar = these.genetic.value.var.covar,
    breeding.value.var.covar = these.breeding.value.var.covar,
    initial.count.value = length(these.genetic.value.var.covar)
  )
   
  #Join the variance-covariance matrices together
  these.trait.var.covar <-  directional.selection.population.10.pct$trait.var.covar
  these.genetic.value.var.covar <- directional.selection.population.10.pct$genetic.value.var.covar
  these.breeding.value.var.covar <- directional.selection.population.10.pct$breeding.value.var.covar
  
  #Join the trait values together
  these.trait.values <- directional.selection.population.10.pct$trait.values
  these.genetic.values <- directional.selection.population.10.pct$genetic.values
  these.breeding.values <- directional.selection.population.10.pct$breeding.values
  
  
  #Simulate this.nGenerations generations of directional selection, selection on 20% of the individuals
  directional.selection.population.20.pct <- cross.stuff.for.a.whole.bunch.of.generations(
    current.generation = TrainPop,
    type.of.selection = "Direct",
    nSelect = 2*this.nSelect,
    nCross = this.nCross,
    nGenerations = this.nGenerations,
    nQtl = this.nQtl,
    SP.within.function = SP,
    trait.var.covar = these.trait.var.covar,
    genetic.value.var.covar = these.genetic.value.var.covar,
    breeding.value.var.covar = these.breeding.value.var.covar,
    initial.count.value = length(these.genetic.value.var.covar)
  )
  
  #Join the variance-covariance matrices together
  these.trait.var.covar <-  directional.selection.population.20.pct$trait.var.covar
  these.genetic.value.var.covar <- directional.selection.population.20.pct$genetic.value.var.covar
  these.breeding.value.var.covar <- directional.selection.population.20.pct$breeding.value.var.covar
  
  #Join the trait values together
  these.trait.values <- directional.selection.population.20.pct$trait.values
  these.genetic.values <- directional.selection.population.20.pct$genetic.values
  these.breeding.values <- directional.selection.population.20.pct$breeding.values
  
  #Simulate this.nGenerations generations of disruptive selection, selection on 10% of the individuals
  disruptive.selection.population.10.pct <- cross.stuff.for.a.whole.bunch.of.generations(
    current.generation = TrainPop,
    type.of.selection = "Disruptive",
    nSelect = this.nSelect,
    nCross = this.nCross,
    nGenerations = this.nGenerations,
    nQtl = this.nQtl,
    SP.within.function = SP,
    trait.var.covar = these.trait.var.covar,
    genetic.value.var.covar = these.genetic.value.var.covar,
    breeding.value.var.covar = these.breeding.value.var.covar,
    initial.count.value = length(these.genetic.value.var.covar)
  )
  
  #Join the variance-covariance matrices together
  these.trait.var.covar <- disruptive.selection.population.10.pct$trait.var.covar
  these.genetic.value.var.covar <- disruptive.selection.population.10.pct$genetic.value.var.covar
  these.breeding.value.var.covar <- disruptive.selection.population.10.pct$breeding.value.var.covar
  
  #Join the trait values together
  these.trait.values <- disruptive.selection.population.10.pct$trait.values
  these.genetic.values <- disruptive.selection.population.10.pct$genetic.values
  these.breeding.values <- disruptive.selection.population.10.pct$breeding.values
  
  
  #Simulate this.nGenerations generations of disruptive selection, selection on 20% of the individuals
  disruptive.selection.population.20.pct <- cross.stuff.for.a.whole.bunch.of.generations(
    current.generation = TrainPop,
    type.of.selection = "Disruptive",
    nSelect = 2*this.nSelect,
    nCross = this.nCross,
    nGenerations = this.nGenerations,
    nQtl = this.nQtl,
    SP.within.function = SP,
    trait.var.covar = these.trait.var.covar,
    genetic.value.var.covar = these.genetic.value.var.covar,
    breeding.value.var.covar = these.breeding.value.var.covar,
    initial.count.value = length(these.genetic.value.var.covar)
  )
  
  #Join the variance-covariance matrices together
  these.trait.var.covar <- disruptive.selection.population.20.pct$trait.var.covar
  these.genetic.value.var.covar <- disruptive.selection.population.20.pct$genetic.value.var.covar
  these.breeding.value.var.covar <- disruptive.selection.population.20.pct$breeding.value.var.covar
  
  #Join the trait values together
  these.trait.values <- disruptive.selection.population.20.pct$trait.values
  these.genetic.values <- disruptive.selection.population.20.pct$genetic.values
  these.breeding.values <- disruptive.selection.population.20.pct$breeding.values
  
  #Simulate this.nGenerations generations of stabilizing selection, selection on 10% of the individuals
  stabilizing.selection.population.10.pct <- cross.stuff.for.a.whole.bunch.of.generations(
    current.generation = TrainPop,
    type.of.selection = "Stabilizing",
    nSelect = this.nSelect,
    nCross = this.nCross,
    nGenerations = this.nGenerations,
    nQtl = this.nQtl,
    SP.within.function = SP,
    trait.var.covar = these.trait.var.covar,
    genetic.value.var.covar = these.genetic.value.var.covar,
    breeding.value.var.covar = these.breeding.value.var.covar,
    initial.count.value = length(these.genetic.value.var.covar)
  )
  #Join the variance-covariance matrices together
  these.trait.var.covar <-  stabilizing.selection.population.10.pct$trait.var.covar
  these.genetic.value.var.covar <- stabilizing.selection.population.10.pct$genetic.value.var.covar
  these.breeding.value.var.covar <- stabilizing.selection.population.10.pct$breeding.value.var.covar

  #Join the trait values together
  these.trait.values <- stabilizing.selection.population.10.pct$trait.values
  these.genetic.values <- stabilizing.selection.population.10.pct$genetic.values
  these.breeding.values <- stabilizing.selection.population.10.pct$breeding.values

  
  
  #Simulate this.nGenerations generations of stabilizing selection, selection on 20% of the individuals
  stabilizing.selection.population.20.pct <- cross.stuff.for.a.whole.bunch.of.generations(
    current.generation = TrainPop,
    type.of.selection = "Stabilizing",
    nSelect = 2*this.nSelect,
    nCross = this.nCross,
    nGenerations = this.nGenerations,
    nQtl = this.nQtl,
    SP.within.function = SP,
    trait.var.covar = these.trait.var.covar,
    genetic.value.var.covar = these.genetic.value.var.covar,
    breeding.value.var.covar = these.breeding.value.var.covar,
    initial.count.value = length(these.genetic.value.var.covar)
  )
  #Join the variance-covariance matrices together
  these.trait.var.covar <-  stabilizing.selection.population.20.pct$trait.var.covar
  names(these.trait.var.covar) <- paste(names(these.trait.var.covar),".Factor.A.", i,
                                        "Factor.B.", j,
                                        "Factor.C.", k,
                                        "Factor.D.", el,
                                        "Rep.", this.rep, sep = "")
  
  these.genetic.value.var.covar <- stabilizing.selection.population.20.pct$genetic.value.var.covar
  names(these.genetic.value.var.covar) <- paste(names(these.genetic.value.var.covar),".Factor.A.", i,
                                        "Factor.B.", j,
                                        "Factor.C.", k,
                                        "Factor.D.", el,
                                        "Rep.", this.rep, sep = "")
  
  these.breeding.value.var.covar <- stabilizing.selection.population.20.pct$breeding.value.var.covar
  names( these.breeding.value.var.covar) <- paste(names(these.breeding.value.var.covar),".Factor.A.", i,
                                                "Factor.B.", j,
                                                "Factor.C.", k,
                                                "Factor.D.", el,
                                                "Rep.", this.rep, sep = "")
  

  #Join the trait values together
  these.trait.values <- stabilizing.selection.population.20.pct$trait.values
  names(these.trait.values) <- paste(names(these.trait.values),".Factor.A.", i,
                                        "Factor.B.", j,
                                        "Factor.C.", k,
                                        "Factor.D.", el,
                                        "Rep.", this.rep, sep = "")
  
  these.genetic.values <- stabilizing.selection.population.20.pct$genetic.values
  names(these.genetic.values) <- paste(names(these.genetic.values),".Factor.A.", i,
                                     "Factor.B.", j,
                                     "Factor.C.", k,
                                     "Factor.D.", el,
                                     "Rep.", this.rep, sep = "")
  
  these.breeding.values <- stabilizing.selection.population.20.pct$breeding.values
  names(these.breeding.values) <- paste(names(these.breeding.values),".Factor.A.", i,
                                       "Factor.B.", j,
                                       "Factor.C.", k,
                                       "Factor.D.", el,
                                       "Rep.", this.rep, sep = "")
  
  
  #Merge the SNP data together from the three different subpopulations
  directional.subpopulation.10.pct.SNPs <- pullSnpGeno(directional.selection.population.10.pct$current.generation, 
                                           simParam = SP)
 
  directional.subpopulation.20.pct.SNPs <- pullSnpGeno(directional.selection.population.20.pct$current.generation, 
                                                       simParam = SP)
   
  disruptive.subpopulation.10.pct.SNPs <- pullSnpGeno(disruptive.selection.population.10.pct$current.generation, 
                                           simParam = SP)
 
  disruptive.subpopulation.20.pct.SNPs <- pullSnpGeno(disruptive.selection.population.20.pct$current.generation, 
                                               simParam = SP)
  
  stabilizing.subpopulation.10.pct.SNPs <- pullSnpGeno(stabilizing.selection.population.10.pct$current.generation, 
                                               simParam = SP)
  
  stabilizing.subpopulation.20.pct.SNPs <- pullSnpGeno(stabilizing.selection.population.20.pct$current.generation, 
                                                simParam = SP)
  
  combined.subpopulation.SNPs <- rbind(directional.subpopulation.10.pct.SNPs,directional.subpopulation.20.pct.SNPs,
                                       disruptive.subpopulation.10.pct.SNPs,disruptive.subpopulation.20.pct.SNPs,
                                       stabilizing.subpopulation.10.pct.SNPs, stabilizing.subpopulation.20.pct.SNPs)
  
  
  #Merge the QTN data together from the three different subpopulations
  directional.subpopulation.QTNs.10.pct <- pullQtlGeno(pop=directional.selection.population.10.pct$current.generation, 
                                                trait = 1, asRaw = FALSE, 
                                                simParam = SP)
 
  directional.subpopulation.QTNs.20.pct <- pullQtlGeno(pop=directional.selection.population.20.pct$current.generation, 
                                                       trait = 1, asRaw = FALSE, 
                                                       simParam = SP) 
  
  
  disruptive.subpopulation.QTNs.10.pct <- pullQtlGeno(pop=disruptive.selection.population.10.pct$current.generation, 
                                               trait = 1, asRaw = FALSE, 
                                               simParam = SP)
 
  disruptive.subpopulation.QTNs.20.pct <- pullQtlGeno(pop=disruptive.selection.population.20.pct$current.generation, 
                                                      trait = 1, asRaw = FALSE, 
                                                      simParam = SP) 
  
  stabilizing.subpopulation.QTNs.10.pct <- pullQtlGeno(pop=stabilizing.selection.population.10.pct$current.generation, 
                                                trait = 1, asRaw = FALSE, 
                                                simParam = SP)
  
  stabilizing.subpopulation.QTNs.20.pct <- pullQtlGeno(pop=stabilizing.selection.population.20.pct$current.generation, 
                                                       trait = 1, asRaw = FALSE, 
                                                       simParam = SP)
  
  combined.subpopulation.QTNs <- rbind(directional.subpopulation.QTNs.10.pct, directional.subpopulation.QTNs.20.pct,
                                       disruptive.subpopulation.QTNs.10.pct,disruptive.subpopulation.QTNs.20.pct,
                                       stabilizing.subpopulation.QTNs.10.pct, stabilizing.subpopulation.QTNs.20.pct)
  
  #Merge the trait data together from the three different subpopulations
  directional.subpopulation.trait.10.pct <- directional.selection.population.10.pct$current.generation@pheno
  row.names(directional.subpopulation.trait.10.pct) = row.names(directional.subpopulation.10.pct.SNPs)

  directional.subpopulation.trait.20.pct <- directional.selection.population.20.pct$current.generation@pheno
  row.names(directional.subpopulation.trait.20.pct) = row.names(directional.subpopulation.20.pct.SNPs)
  
  disruptive.subpopulation.trait.10.pct <- disruptive.selection.population.10.pct$current.generation@pheno
  row.names(disruptive.subpopulation.trait.10.pct) = row.names(disruptive.subpopulation.10.pct.SNPs)

  disruptive.subpopulation.trait.20.pct <- disruptive.selection.population.20.pct$current.generation@pheno
  row.names(disruptive.subpopulation.trait.20.pct) = row.names(disruptive.subpopulation.20.pct.SNPs)  
  
  stabilizing.subpopulation.trait.10.pct <- stabilizing.selection.population.10.pct$current.generation@pheno
  row.names(stabilizing.subpopulation.trait.10.pct) = row.names(stabilizing.subpopulation.10.pct.SNPs)

  stabilizing.subpopulation.trait.20.pct <- stabilizing.selection.population.20.pct$current.generation@pheno
  row.names(stabilizing.subpopulation.trait.20.pct) = row.names(stabilizing.subpopulation.20.pct.SNPs)  
  
  combined.subpopulation.trait <- rbind(directional.subpopulation.trait.10.pct, directional.subpopulation.trait.20.pct,
                                        disruptive.subpopulation.trait.10.pct,disruptive.subpopulation.trait.20.pct,
                                        stabilizing.subpopulation.trait.10.pct,stabilizing.subpopulation.trait.20.pct)
  
  
  ###Save the R workspace
  #I am not sure if the code starting on Line 273 will work. Thus,
  # I want to save the image here. Once I am convinced this is 
  # correctly, I will comment this out
  save.image(paste(factor.A,".FactorA..",factor.B, ".FactorB..",
                   factor.C,".FactorC..",factor.D,".FactorD..",
                   this.rep, ".Rep.Rdata", sep = ""))
  
  
  #Merge the trait data together from the three different subpopulations
  directional.subpopulation.trait.10.prev.gen <- directional.selection.population.10.pct$previous.generation@pheno
  row.names(directional.subpopulation.trait.10.prev.gen) = row.names(pullSnpGeno(directional.selection.population.10.pct$previous.generation, simParam = SP))
  

  directional.subpopulation.trait.20.prev.gen <- directional.selection.population.20.pct$previous.generation@pheno
  row.names(directional.subpopulation.trait.20.prev.gen) = row.names(pullSnpGeno(directional.selection.population.20.pct$previous.generation, simParam = SP))
  
  disruptive.subpopulation.trait.10.prev.gen <- disruptive.selection.population.10.pct$previous.generation@pheno
  row.names(disruptive.subpopulation.trait.10.prev.gen) = row.names(pullSnpGeno(disruptive.selection.population.10.pct$previous.generation, simParam = SP))
  
  disruptive.subpopulation.trait.20.prev.gen <- disruptive.selection.population.20.pct$previous.generation@pheno
  row.names(disruptive.subpopulation.trait.20.prev.gen) = row.names(pullSnpGeno(disruptive.selection.population.20.pct$previous.generation, simParam = SP))  
  
  stabilizing.subpopulation.trait.10.prev.gen <- stabilizing.selection.population.10.pct$previous.generation@pheno
  row.names(stabilizing.subpopulation.trait.10.prev.gen) = row.names(pullSnpGeno(stabilizing.selection.population.10.pct$previous.generation, simParam = SP))
  
  stabilizing.subpopulation.trait.20.prev.gen <- stabilizing.selection.population.20.pct$previous.generation@pheno
  row.names(stabilizing.subpopulation.trait.20.prev.gen) = row.names(pullSnpGeno(stabilizing.selection.population.20.pct$previous.generation, simParam = SP))  
  
  combined.subpopulation.trait <- rbind(directional.subpopulation.trait.10.prev.gen, directional.subpopulation.trait.20.prev.gen,
                                        disruptive.subpopulation.trait.10.prev.gen,disruptive.subpopulation.trait.20.prev.gen,
                                        stabilizing.subpopulation.trait.10.prev.gen,stabilizing.subpopulation.trait.20.prev.gen)
  
  save.image(paste(factor.A,".FactorA..",factor.B, ".FactorB..",
                   factor.C,".FactorC..",factor.D,".FactorD..",
                   this.rep, ".Rep.Rdata", sep = ""))
  
#}#End one.replicate.of.simulation.run
