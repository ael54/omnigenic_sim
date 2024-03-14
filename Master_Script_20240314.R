#Please see Slide 8: of /Lab_Presentations/Lab_Meeting_20240208/Lipka_A_Omnigenic_20240208.pptx
# for a visualization of what I want to do with this simualtion

#Clear your workspace
rm(list = ls())

#Set the working directory
setwd("/Users/alipka/Library/CloudStorage/Box-Box/Sabbatical_Roslin_Institute/R_workspace/Sabbatical_Project")
home.dir <- getwd()

#Open all libraries
library(package = "AlphaSimR")
library(package = "simplePHENOTYPES")

#Source in all files
source("Simulating_Omnigenic_Genetic_Architecture/Simulate_Omnigenic_Genetic_Architecture_as_a_Function_20240219.R")
source("Simulating_Omnigenic_Genetic_Architecture/Simulate_Multiple_Generations_of_Selection.R")
source("Functions_to_Make_Life_Easier/get.me.my.SNPs.in.hapmap.format.R")
################################################################################
################################################################################
#Simulate a founder population
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

####This creates the founder population, which is necessary for simulating traits.
#the.founders <- newPop(founderPop)


#Simulate an initial genetic architecture in the founder population

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


stabilizing.subpopulation.QTNs <- pullQtlGeno(pop=disruptive.selection.population, 
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
