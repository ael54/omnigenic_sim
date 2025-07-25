#Set the working directory

#The first time you run this, please be sure to install
# AlphaSimR and simplePHENOTYPES



#Set your working directory
setwd("/Users/alipka/Library/CloudStorage/Box-Box/Sabbatical_Roslin_Institute/R_workspace/omnigenic_sim/")
home.dir <- getwd()


#Open all libraries
library(package = "AlphaSimR")
library(package = "simplePHENOTYPES")

#Source in all files

source("Simulating_Omnigenic_Genetic_Architecture/Simulate_Omni_Four_Trait_Trick_20240509.R")
source("Simulating_Omnigenic_Genetic_Architecture/Simulate_Multiple_Generations_of_Selection_20241216.R")
source("Simulating_Omnigenic_Genetic_Architecture/Obtain_Breeding_Values_20240515.R")
source("Simulating_Omnigenic_Genetic_Architecture/Simulate_Omnigenic_Genetic_Architecture_as_a_Function_20240219.R")
source("Functions_to_Make_Life_Easier/get.me.my.SNPs.in.hapmap.format.R")
source("Simulate_a_founder_Population_as_a_Function_20250610.R")

#Simulate a founder population in AlphaSimR. This will take a long-ish time to run
# e.g. - it should take about 30 minutes on my MacBook Pro


###For loop through all settings of the factorial experiment
### The end product should be a whole bunch of R workspaces
### Not to mention output directories of simulated traits

for(rep in 1:3){
  simulated.founder.pop.genos <- simulate.a.founder.population(this.nChr = 2,
                                                               this.nQtl = 50,
                                                               this.nSnp = 400)
  this.QTL.Map <- simulated.founder.pop.geno$this.QTL.Map
  these.QTL.Genotypes <- simulated.founder.pop.geno$these.QTL.Genotypes
  the.physical.map.of.QTLs <- simulated.founder.pop.geno$the.physical.map.of.QTLs
  hapmap.file.of.founder.QTLs <- simulated.founder.pop.geno$hapmap.file.of.founder.QTLs
  the.founder.SNPs <- simulated.founder.pop.geno$the.founder.SNPs
  the.physical.map.of.SNPs <- simulated.founder.pop.geno$the.physical.map.of.SNPs
  hapmap.file.of.founder.SNPs <- simulated.founder.pop.geno$hapmap.file.of.founder.SNPs
  
  for(i in c(1,2,4)){
    for(j in c(0.05,0.5,1,2)){
       for(k in c(0.05,0.5,1,2)){
         for(el in c(0.05,0.5,1,2)){
           #Run the function
           factor.A <- i
           factor.B <- j
           factor.C <- k
           factor.D <- el 
           this.rep <- rep
           this.output.directory.name <-  paste("Factor.A.", i,
                                                "Factor.B.", j,
                                                "Factor.C.", k,
                                                "Factor.D.", el,
                                                "Rep.", this.rep,
                                                sep = "")
           source("Master_Simulation_Script_20241216.R")
         }#End for(rep in 1:3)
      }#End for(el in 1:4)
    }#End for(k in 1:4)
  }#End for(j in 1:4)
}#End for(i in 1:3)

