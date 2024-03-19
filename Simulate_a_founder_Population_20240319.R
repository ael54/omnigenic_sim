#Please see Slide 8: of /Lab_Presentations/Lab_Meeting_20240208/Lipka_A_Omnigenic_20240208.pptx
# for a visualization of what I want to do with this simualtion




  #Source in all files
  source("Simulating_Omnigenic_Genetic_Architecture/Simulate_Omnigenic_Genetic_Architecture_as_a_Function_20240219.R")
  source("Simulating_Omnigenic_Genetic_Architecture/Simulate_Multiple_Generations_of_Selection.R")
  source("Functions_to_Make_Life_Easier/get.me.my.SNPs.in.hapmap.format.R")
  ################################################################################
  ################################################################################
  #Simulate a founder population
  print("-------------------All simulation settings for founder popuations are found on lines 16-55 of Simulate_a_founder_Population_20240319--------------------")
  
  this.nChr <- 5
  this.nQtl <- 50  # This is actually the number of QTLs per chromosome
  this.nSnp <- 1000 # This is the total number of SNPs across the genome
  this.Ne <- 30
  this.nParents <- 50
  #this.nInd <- 300 #This is obsolete because we are using this.nParents instead
  this.histNe.vector <- c(100, 100000, 1000000) 
  this.histGen.vector <- c(100, 1000, 10000)
  this.split <- NULL
  this.nSelect <- 112 # This is based on selecting 10% of individuals in the training population 
  this.nCross <- 1120 #200 is probably too small to be realistic
  this.nGenerations <- 10 #5 is probably too small to be realistic
  ######
  ##Input parameters for the burnin phase. These numbers are based on Bancic et al. (2023)
  this.varE.during.burnin <- 4
  this.repEYT <- 8
  # Other input parameters for the burnin phase, to be used during the 
  #.  fill in stage
  # ---- Breeding program details for developing the founder population----
  nParents = this.nParents  # Number of parents to start a breeding cycle
  nCrosses = this.nCross # Number of crosses per year
  nDH      = 100 # DH lines produced per cross
  famMax   = 10  # The maximum number of DH lines per cross to enter PYT
  nPYT     = 500 # Entries per preliminary yield trial
  nAYT     = 50  # Entries per advanced yield trial
  nEYT     = 10  # Entries per elite yield trial
  varE     = 4
  
  # Effective replication of yield trials
  repHDRW  = 4/9 # h2 = 0.1
  repPYT   = 1   # h2 = 0.2
  repAYT   = 4   # h2 = 0.5
  repEYT   = 8   # h2 = 0.7
  
  #Add information needed to do the Burn-in
  nReps   = 1    # Number of simulation replicates
  nBurnin = 20   # Number of years in burnin phase
  nFuture = 20   # Number of years in future phase
  nCycles = nBurnin + nFuture
  startTP = 19   # Year to start training population
  
  
  
  #The comments below are paraphrased from the AlphaSim demos
  
  ###############################################################################
  ###############################################################################
  ###############################################################################
  #Create a founder population
  ############## Modification on March 6, 2024. I am going to do 
  ### a fill-in and burn-in phase so that we can have genomic
  ### properties (e.g., LD) that resemble what we expect to see
  ### in a founder population. I am going to use the parameters
  ### that were set in Bancic et al. (2023). The rationale is that
  ### I (Alex Lipka) am not an expert in this, and I am sure that
  #### these starting parameters were pressure tested.
  #### All of these are based on/use code from "
  ###### jbancic_alphasimr_plants/01_LineBreeding/03_TwoPartGS/"
  #This code is a slightly modified version of CreateParents.R 
  #For the time being, I am only going to simulate one rep. 
  # If we end up wanting to simulate more than one rep, a for loop
  # will be put around the code below, as done beginning on 
  # Line 25 of 00RUNME.R. If we decide to do this, replace
  # "REP <- 1" below with the beginning of the for loop
  REP <- 1
  output = data.frame(year     = 1:nCycles,
                      rep      = rep(REP, nCycles),
                      scenario = "",
                      meanG = numeric(nCycles),
                      varG  = numeric(nCycles),
                      accSel  = numeric(nCycles))
  
  #founderPop = runMacs2(nInd = this.nParents,
   #                    nChr = this.nChr,
  #                     segSites = this.nQtl +this.nSnp,
  #                     Ne = this.Ne,
  #                    histNe = this.histNe.vector,
  #                     histGen = this.histGen.vector,
  #                     split = this.split,
  #                     inbred = TRUE)
  
  founderPop = quickHaplo(nInd = this.nParents,
                        nChr = this.nChr,
                        segSites = this.nQtl +this.nSnp,
                        #Ne = this.Ne,
                        #histNe = this.histNe.vector,
                        #histGen = this.histGen.vector,
                        #split = this.split,
                        inbred = TRUE)
  
  #Set population parameters
  SP = SimParam$new(founderPop)
  
  # Add SNP chip
  SP$restrSegSites(this.nQtl, this.nSnp)
  if (this.nSnp > 0) {
    SP$addSnpChip(nSnpPerChr = this.nSnp/this.nChr)
  }
  
  #The file sourced in below is a slightly modified version of "CreateParents.R"
  source(file = "Slightly_tweaked_scripts_from_Bancic_et_al/Finish_Creating_the_Parents.R") 
  
  
  #Thie file sourced in below is a slightly modified version of FillPipeline.R 
  
  source(file = "Slightly_tweaked_scripts_from_Bancic_et_al/FillPipeline_Minimal_Mods.R") 
  
  
  #Now, do the burn-in period. This code is from Lines 46-55 of 00RUNME.R
  cat("--> Working on Burn-in \n")
  for(year in 1:nBurnin) {
    cat(" Working on burn-in year:",year,"\n")
    source(file = "Slightly_tweaked_scripts_from_Bancic_et_al/UpdateParents_mod.R") # Pick new parents
    source(file = "Slightly_tweaked_scripts_from_Bancic_et_al/AdvanceYear_mod.R")   # Advance yield trials by a year
    source(file = "Slightly_tweaked_scripts_from_Bancic_et_al/StoreTrainPop_mod.R") # Store training population
    # Report results
    output$meanG[year] = meanG(DH)
    output$varG[year]  = varG(DH)
  }
  
  #At the end of this code, the founder population will be called "TrainPop"
  # This is an "S4 object of class Pop"
  
  ###########################################################################
  ###########################################################################
  ###########################################################################
  ###########################################################################
  ###########################################################################
  ###########################################################################
  
  
  #Now, we start simulating a trait following a putative omnigenic genetic
  # architecture in this founder population. 
  
  
  #Get the QTLs
  this.QTL.Map <- getQtlMap(trait = 1, simParam = SP)
  these.QTL.Genotypes <- pullQtlGeno(pop=TrainPop, trait = 1, asRaw = FALSE, 
                                     simParam = SP)
  
  
  the.physical.map.of.QTLs <- getQtlMap(trait = 1, simParam = SP)
  
  #Format the QTLs for reading into simplePHENOTYPES...aw shucks!
  hapmap.file.of.founder.QTLs <- get.me.my.SNPs.in.hapmap.format(these.SNPs = these.QTL.Genotypes,
                                                                 this.physical.map = this.QTL.Map)
  
  
  
  #Get your SNPs
  the.founder.SNPs <- pullSnpGeno(TrainPop, simParam = SP)
  
  
  #Get your genetic map of the SNPs
  the.physical.map.of.SNPs <- getSnpMap(snpChip = 1, simParam = SP)
  
  #Prepare everything for reading into simplePHENOTYPES
  hapmap.file.of.founder.SNPs <- get.me.my.SNPs.in.hapmap.format(these.SNPs = the.founder.SNPs,
                                  this.physical.map = the.physical.map.of.SNPs)
  
  

  