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
this.nParents <- 10
this.nInd <- 300 #300 is probably too small to be realistic
this.histNe.vector <- c(100, 100000, 1000000) 
this.histGen.vector <- c(100, 1000, 10000)
this.split <- NULL
this.nSelect <- 25 #25 is probably too small to be realistic 
this.nCross <- 200 #200 is probably too small to be realistic
this.nGenerations <- 1 #5 is probably too small to be realistic
######
##Input parameters for the burnin phase
this.varE.during.burnin <- 4
this.repEYT <- 8
# Other input parameters for the burnin phase, to be used during the 
#.  fill in stage
# ---- Breeding program details ----
nParents = 50  # Number of parents to start a breeding cycle
nCrosses = 100 # Number of crosses per year
nDH      = 100 # DH lines produced per cross
famMax   = 10  # The maximum number of DH lines per cross to enter PYT
nPYT     = 500 # Entries per preliminary yield trial
nAYT     = 50  # Entries per advanced yield trial
nEYT     = 10  # Entries per elite yield trial


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
founderPop = runMacs2(nInd = this.nParents,
                     nChr = this.nChr,
                     segSites = this.nQtl +this.nSnp,
                     Ne = this.Ne,
                     histNe = this.histNe.vector,
                     histGen = this.histGen.vector,
                     split = this.split,
                     inbred = TRUE)

#Set population parameters
SP = SimParam$new(founderPop)

# Add SNP chip
SP$restrSegSites(this.nQtl, this.nSnp)
if (this.nSnp > 0) {
  SP$addSnpChip(nSnpPerChr = this.nSnp/this.nChr)
}


#Simulate a trait. We are not actually going to use this trait for anything
# once we have a founder population, so I will hard code in the settings
# from Bancic et al. (2023)
SP$addTraitAG(nQtlPerChr = this.nQtl,
              mean       = 1,
              var        = 1,
              varEnv     = 1,
              varGxE     = 4)

# Collect pedigree
SP$setTrackPed(TRUE)

###This code below is from the wheat example of Bancic et al. (2023)
# Create founder parents
Parents = newPop(founderPop)

# Add phenotype reflecting evaluation in EYT
Parents = setPheno(Parents, varE = this.varE.during.burnin, reps = this.repEYT)






#Add a SNP chip


#Add a dummy trait that is not actually used. This will let us pass the 
# QTLs to simplePHENOTYPES, and we can run downstream analysis on SNPs
# instead of QTLs





######### Do this later after you did the burn ins
#Get the QTLs
this.QTL.Map <- getQtlMap(trait = 1, simParam = SP)
these.QTL.Genotypes <- pullQtlGeno(pop=founderPop, trait = 1, asRaw = FALSE, 
                                   simParam = SP)


the.physical.map.of.QTLs <- getQtlMap(trait = 1, simParam = SP)

#Format the QTLs for reading into simplePHENOTYPES...aw shucks!
hapmap.file.of.founder.QTLs <- get.me.my.SNPs.in.hapmap.format(these.SNPs = these.QTL.Genotypes,
                                                               this.physical.map = this.QTL.Map)



#Get your SNPs
the.founder.SNPs <- pullSnpGeno(the.founders, simParam = SP)


#Get your genetic map of the SNPs
the.physical.map.of.SNPs <- getSnpMap(snpChip = 1, simParam = SP)

#Prepare everything for reading into simplePHENOTYPES
hapmap.file.of.founder.SNPs <- get.me.my.SNPs.in.hapmap.format(these.SNPs = the.founder.SNPs,
                                this.physical.map = the.physical.map.of.SNPs)

####This creates the founder population, which is necessary for simulating traits.
the.founders <- newPop(founderPop)


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
the.founders@pheno <- as.matrix(this.simulated.trait$this.simulated.trait[,2])
  
################################################################################
################################################################################
#For loop through various breeding programs (directional seletion,
#            disruptive selection, stabilizing selection)... 

  
#Simulate this.nGenerations generations of directional selection
directional.selection.population <- cross.stuff.for.a.whole.bunch.of.generations(
  current.generation = the.founders,
  type.of.selection = "Direct",
  nSelect = this.nSelect,
  nCross = this.nCross,
  nGenerations = this.nGenerations,
  nQtl = this.nQtl,
  SP.within.function = SP
)
 

#Simulate this.nGenerations generations of disruptive selection
disruptive.selection.population <- cross.stuff.for.a.whole.bunch.of.generations(
  current.generation = the.founders,
  type.of.selection = "Disruptive",
  nSelect = this.nSelect,
  nCross = this.nCross,
  nGenerations = this.nGenerations,
  nQtl = this.nQtl,
  SP.within.function = SP
)


#Simulate this.nGenerations generations of stabilizing selection
stabilizing.selection.population <- cross.stuff.for.a.whole.bunch.of.generations(
  current.generation = the.founders,
  type.of.selection = "Stabilizing",
  nSelect = this.nSelect,
  nCross = this.nCross,
  nGenerations = this.nGenerations,
  nQtl = this.nQtl,
  SP.within.function = SP
)



#Merge the genotype data together from the three different subpopulations
directional.subpopulation.SNPs <- pullSnpGeno(directional.selection.population, 
                                         simParam = SP)

disruptive.subpopulation.SNPs <- pullSnpGeno(disruptive.selection.population, 
                                         simParam = SP)

stabilizing.subpopulation.SNPs <- pullSnpGeno(stabilizing.selection.population, 
                                             simParam = SP)

combined.subpopulation.SNPs <- rbind(directional.subpopulation.SNPs,
                          disruptive.subpopulation.SNPs,stabilizing.subpopulation.SNPs)


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
save.image("Test_Run_small_data_set_20240305.Rdata")
