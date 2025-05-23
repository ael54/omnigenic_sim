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
this.nChr <- 2
this.nQtl <- 0 #I am not going to have QTLs here because simplePHENOTYPES will 
               # the SNPs as QTLs (which we call QTNs)
this.nSnp <- 250
this.nInd <- 300
this.nSelect <- 25 
this.nCross <- 200 
this.nGenerations <- 5

#The comments below are paraphrased from the AlphaSim demos

#Create a founder population
founderPop = runMacs(nInd = this.nInd,
                     nChr = this.nChr,
                     segSites = this.nQtl +this.nSnp,
                     inbred = TRUE,
                     species = "MAIZE")

#Set population parameters
SP = SimParam$new(founderPop)


#Add a SNP chip
SP$addSnpChip(nSnpPerChr = this.nSnp/this.nChr)


#Obtain the founder population
the.founders <- newPop(founderPop)

#Get your SNPs
the.founder.SNPs <- pullSnpGeno(the.founders, simParam = SP)


#Get your genetic map of the SNPs
the.physical.map.of.SNPs <- getSnpMap(snpChip = 1, simParam = SP)

#Prepare everything for reading into simplePHENOTYPES
hapmap.file.of.founder.SNPs <- get.me.my.SNPs.in.hapmap.format(these.SNPs = the.founder.SNPs,
                                this.physical.map = the.physical.map.of.SNPs)


 

################################################################################
################################################################################
#Simulate an initial genetic architecture in the founder population

#Initialize the parameters for simulating the omnigenic genetic architecture the function

source("Simulating_Omnigenic_Genetic_Architecture/setInputParametersforSimplePhentoypes.R")
this.input.SNPs <- hapmap.file.of.founder.SNPs


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
  SP.within.function = SP
)
 

#Simulate this.nGenerations generations of disruptive selection
disruptive.selection.population <- cross.stuff.for.a.whole.bunch.of.generations(
  current.generation = the.founders,
  type.of.selection = "Disruptive",
  nSelect = this.nSelect,
  nCross = this.nCross,
  nGenerations = this.nGenerations,
  SP.within.function = SP
)


#Simulate this.nGenerations generations of stabilizing selection
stabilizing.selection.population <- cross.stuff.for.a.whole.bunch.of.generations(
  current.generation = the.founders,
  type.of.selection = "Stabilizing",
  nSelect = this.nSelect,
  nCross = this.nCross,
  nGenerations = this.nGenerations,
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
directional.subpopulation.trait <- directional.selection.population@pheno
row.names(directional.subpopulation.trait) = row.names(directional.subpopulation.SNPs)

disruptive.subpopulation.trait <- disruptive.selection.population@pheno
row.names(disruptive.subpopulation.trait) = row.names(disruptive.subpopulation.SNPs)


stabilizing.subpopulation.trait <- stabilizing.selection.population@pheno
row.names(stabilizing.subpopulation.trait) = row.names(stabilizing.subpopulation.SNPs)


combined.subpopulation.trait <- rbind(directional.subpopulation.trait,
                  disruptive.subpopulation.trait,stabilizing.subpopulation.trait)


###Save the R workspace
save.image("Test_Run_small_data_set_20240223.Rdata")
