#Set your working directory
setwd("/Users/alipka/Library/CloudStorage/Box-Box/Sabbatical_Roslin_Institute/R_workspace/omnigenic_sim")
home.dir <- getwd()

#####Read in all of the packages that are necessary
#Read in prerequiste libaries for GAPIT
library('MASS')
library(multtest)
library(gplots)
library(sommer)
library(AlphaSimR)



#Read in GAPIT
setwd("Scripts_Necessary_for_GAPIT")
source("GAPIT_EMMA source code.txt")
source("GAPIT_Code_from_Internet_20120411_Allelic_Effect.r")
setwd(home.dir)



#This count variable will help with appending to these lists storing
# Multi_Kernel_BLUP GS summary statistics
count <- 1

#Initiallize objects that will store GS results
these.prediction.accuracies.SNPs.Multi.Kernel.BLUP__ADD_only <- list()

###For loop through all settings of the factorial experiment
### The end products populated versions of the updated lists
#### where each item of each list corresponds to results
#### from a particular setting of the factorial experiment
### Not to mention output directories of simulated traits

#Alex added this new line below to point to the directory where Joao's run of
# one iteration of the loop was kept
setwd("/Users/alipka/Library/CloudStorage/Box-Box/IR-281/Partial-Results-Simulations-20240918/")

for(i in c(1,2,4)){
  for(j in c(0.05,0.5,1,2)){
    for(k in c(0.05,0.5,1,2)){
       for(el in c(0.05,0.5,1,2)){
         for(rep in 1:3){
           #Run the function
           factor.A <- i
           factor.B <- j
           factor.C <- k
           factor.D <- el 
           this.rep <- rep
           #Read in the R data file for a given setting
           load(paste(factor.A,".FactorA..",factor.B, ".FactorB..",
                 factor.C,".FactorC..",factor.D,".FactorD..",
                 this.rep, ".Rep.Rdata", sep = ""))
           #Keep track of what setting you are on
           this.setting <-  paste(factor.A,".FactorA..",factor.B, ".FactorB..",
                                                 factor.C,".FactorC..",factor.D,".FactorD..",
                                                 this.rep, ".Rep", sep = "")

           #Run the master script; source it in
           source("Master_Multi_Kernel_BLUP_Analysis_Script_ADD_Only_SNPs.R")
           
           ##############################################################
           #### Append GS Results
           #Append results to the lists you created before the start of this rep
           these.prediction.accuracies.SNPs.Multi.Kernel.BLUP__ADD_only[[count]] <- prediction.accuracies.SNPs
           names(these.prediction.accuracies.SNPs.Multi.Kernel.BLUP__ADD_only)[count] <- this.setting 
           #Update count for the next iteration of this loop
           count <- count + 1
         }#End for(rep in 1:3)
      }#End for(el in 1:4)
    }#End for(k in 1:4)
  }#End for(j in 1:4)
}#End for(i in 1:3)

