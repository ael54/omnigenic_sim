#Please see Slide 8: of /Lab_Presentations/Lab_Meeting_20240208/Lipka_A_Omnigenic_20240208.pptx
# for a visualization of what I want to do with this simualtion


simulate.a.founder.population <- function(this.nChr = 2, this.nQtl = 50,  this.nSnp = 400,
                                          this.Ne = 30,
                                          this.nParents = 50,
                                          this.histNe.vector = c(100, 100000, 1000000),
                                          this.histGen.vector = c(100, 1000, 10000),
                                          this.split = NULL,
                                          this.nSelect = 112,  
                                          this.nCross = 1120,
                                          this.nGenerations = 10,
                                          this.varE.during.burnin = 4,
                                          nDH = 100,
                                          famMax = 10,
                                          nPYT = 500,
                                          nAYT = 50, 
                                          nEYT = 10,
                                          varE = 4,
                                          repHDRW  = 4/9, 
                                          repPYT = 1,
                                          repAYT = 4,
                                          repEYT = 8,
                                          nReps   = 1,
                                          nBurnin = 20,
                                          nFuture = 2,
                                          startTP = 19) 
{    
    nCycles = nBurnin + nFuture
    nParents = this.nParents  # Number of parents to start a breeding cycle
    nCrosses = this.nCross # Number of crosses per year

print(paste("famMax = ", famMax, sep = ""))
print(paste("this.nSelect = ", this.nSelect, sep = ""))
print(paste("nCycles = ", nCycles, sep = ""))
  #Source in all files

  ################################################################################
  ################################################################################
  #Simulate a founder population
  
  
  
  
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
  
  print("Hey this is just before REP")
  REP <- 1
  print("Hey this is just before data frame ")
  output = data.frame(year     = 1:nCycles,
                      rep      = rep(REP, nCycles),
                      scenario = "",
                      meanG = numeric(nCycles),
                      varG  = numeric(nCycles),
                      accSel  = numeric(nCycles))
  print("Hey it made it to past data frame ")
  
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
  print("Hey it made it to founder pop")
  #Set population parameters
  SP = SimParam$new(founderPop)
  
  # Add SNP chip
  SP$restrSegSites(this.nQtl, this.nSnp)
  if (this.nSnp > 0) {
    SP$addSnpChip(nSnpPerChr = this.nSnp/this.nChr)
  }
  
  print("Hey it made it through all the stuff before being sourced")
  #The file sourced in below is a slightly modified version of "CreateParents.R"
  #source(file = "Slightly_tweaked_scripts_from_Bancic_et_al/Finish_Creating_the_Parents.R") 
  #Simulate a trait. We are not actually going to use this trait for anything
  # once we have a founder population, so I will hard code in the settings
  # from Bancic et al. (2023)
  SP$addTraitAG(nQtlPerChr = this.nQtl,
                mean       = 1,
                var        = 1,
                varEnv     = 1,
                varGxE     = 4)
  
  print("Oh my gosh. It made it to Line 119")
  # Collect pedigree
  SP$setTrackPed(TRUE)
  print("Oh my gosh. It made it to Line 121") 
  
  # Create founder parents
  Parents = newPop(founderPop)
  
  print("Oh my gosh. It made it to Line 127") 
  # Add phenotype reflecting evaluation in EYT
  Parents = setPheno(Parents, varE = this.varE.during.burnin, reps = this.repEYT)

  print("Oh my gosh. It made it to Line 130") 
  #Thie file sourced in below is a slightly modified version of FillPipeline.R 
  
  #source(file = "Slightly_tweaked_scripts_from_Bancic_et_al/FillPipeline_Minimal_Mods.R") 
  
  for(cohort in 1:7){
    cat("  FillPipeline stage:",cohort,"of 7\n")
    if(cohort<7){
      # Stage 1
      F1 = randCross(Parents, nCrosses)
    }
    if(cohort<6){
      # Stage 2
      DH = makeDH(F1, nDH)
    }
    if(cohort<5){
      # Stage 3
      HDRW = setPheno(DH, varE = varE, reps = repHDRW)
    }
    if(cohort<4){
      # Stage 4
      PYT = selectWithinFam(HDRW, famMax)
      PYT = selectInd(PYT, nPYT)
      PYT = setPheno(PYT, varE = varE, reps = repPYT)
    }
    if(cohort<3){
      # Stage 5
      AYT = selectInd(PYT, nAYT)
      AYT = setPheno(AYT, varE = varE, reps = repAYT)
    }
    if(cohort<2){
      # Stage 6
      EYT = selectInd(AYT, nEYT)
      EYT = setPheno(EYT, varE = varE, reps = repEYT)
    }
    if(cohort<1){
      # Stage 7
    }
  }
  
  #Now, do the burn-in period. This code is from Lines 46-55 of 00RUNME.R
  cat("--> Working on Burn-in \n")
  for(year in 1:nBurnin) {
    cat(" Working on burn-in year:",year,"\n")
    #source(file = "Slightly_tweaked_scripts_from_Bancic_et_al/UpdateParents_mod.R") # Pick new parents
    # Update parents
    
    # Replace 10 oldest inbred parents with 10 new parents from EYT stage
    Parents = c(Parents[11:nParents], EYT)
    
    #source(file = "Slightly_tweaked_scripts_from_Bancic_et_al/AdvanceYear_mod.R")   # Advance yield trials by a year
    # Advance year
    
    # Advance breeding program by 1 year
    # Works backwards through pipeline to avoid copying data
    
    # Stage 7
    # Release variety
    
    # Stage 6
    EYT = selectInd(AYT, nEYT)
    EYT = setPheno(EYT, varE = varE, reps = repEYT)
    
    # Stage 5
    AYT = selectInd(PYT, nAYT)
    AYT = setPheno(AYT, varE = varE, reps = repAYT)
    
    # Stage 4
    output$accSel[year] = cor(HDRW@gv, HDRW@pheno)
    PYT = selectWithinFam(HDRW, famMax)
    PYT = selectInd(PYT, nPYT)
    PYT = setPheno(PYT, varE = varE, reps = repPYT)
    
    # Stage 3
    HDRW = setPheno(DH, varE = varE, reps = repHDRW)
    
    # Stage 2
    DH = makeDH(F1, nDH)
    
    # Stage 1
    F1 = randCross(Parents, nCrosses)
    
    #source(file = "Slightly_tweaked_scripts_from_Bancic_et_al/StoreTrainPop_mod.R") # Store training population
    # Report results
    # Store training population
    
    if (year == startTP){
      cat("  Start collecting training population \n")
      TrainPop = c(PYT, EYT, AYT)
    }
    
    if (year > startTP & year < nBurnin+1){
      cat("  Collecting training population \n")
      TrainPop = c(TrainPop,
                   PYT, EYT, AYT)
    }
    
    if (year > nBurnin){
      cat("  Maintaining training population \n")
      TrainPop = c(TrainPop[-c(1:c(PYT, EYT, AYT)@nInd)],
                   PYT, EYT, AYT)
    }
    
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
  
  return(list(this.QTL.Map = this.QTL.Map, these.QTL.Genotypes = these.QTL.Genotypes,
              the.physical.map.of.QTLs = the.physical.map.of.QTLs, 
              hapmap.file.of.founder.QTLs = hapmap.file.of.founder.QTLs,
              the.founder.SNPs = the.founder.SNPs, the.physical.map.of.SNPs = the.physical.map.of.SNPs,
              hapmap.file.of.founder.SNPs =  hapmap.file.of.founder.SNPs))
}# end simulate.a.founder.population
  