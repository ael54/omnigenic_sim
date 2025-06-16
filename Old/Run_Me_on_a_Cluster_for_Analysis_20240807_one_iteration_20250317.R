


#Set your working directory
setwd("/Users/alipka/Library/CloudStorage/Box-Box/Sabbatical_Roslin_Institute/R_workspace/omnigenic_sim")
home.dir <- getwd()

#####Read in all of the packages that are necessary
#Read in prerequiste libaries for GAPIT
library('MASS')
library(multtest)
library(gplots)
library(sommer)



#Read in GAPIT
setwd("Scripts_Necessary_for_GAPIT")
source("GAPIT_EMMA source code.txt")
source("GAPIT_Code_from_Internet_20120411_Allelic_Effect.r")
setwd(home.dir)



#This count variable will help with appending to these lists storing
# GWAS and GS summary statistics
count <- 1

#Initialize objects that will store GWAS results
these.spearman.rank.correlation.between.GWAS.core.QTNs <- list()
these.median.spearman.correlations.between.core.QTN <- list()
these.median.core.spearman.correlation.within.types.of.selection.QTN <- list()
these.median.core.spearman.correlation.different.types.of.selection.QTN <- list()

these.spearman.rank.correlation.between.GWAS.peripheral.QTNs <- list()
these.median.spearman.correlations.between.peripheral.QTN <- list()
these.median.peripheral.spearman.correlation.within.types.of.selection.QTN <- list()
these.median.peripheral.spearman.correlation.different.types.of.selection.QTN <- list()

these.spearman.rank.correlation.between.GWAS.core.SNPs <- list()
these.median.spearman.correlations.between.core.SNPs <- list()
these.median.core.spearman.correlation.within.types.of.selection.SNPs <- list()
these.median.core.spearman.correlation.different.types.of.selection.SNPs <- list()

these.spearman.rank.correlation.between.GWAS.peripheral.SNPs <- list()
these.median.spearman.correlations.between.peripheral.SNPs <- list()
these.median.peripheral.spearman.correlation.within.types.of.selection.SNPs <- list()
these.median.peripheral.spearman.correlation.different.types.of.selection.SNPs <- list()

#Initiallize objects that will store GS results
these.prediction.accuracies.QTNs <- list()
these.prediction.accuracies.SNPs <- list()

###For loop through all settings of the factorial experiment
### The end products populated versions of the updated lists
#### where each item of each list corresponds to results
#### from a particular setting of the factorial experiment
### Not to mention output directories of simulated traits

#Alex added this new line below to point to the directory where Joao's run of
# one iteration of the loop was kept
setwd("/Users/alipka/Library/CloudStorage/Box-Box/IR-281/simulation-results-for-evaluation-20250317/")

for(i in c(1,2,4)){
  for(j in c(0.05,0.5,1,2)){
    for(k in c(0.05,0.5,1,2)){
       for(el in c(0.05,0.5,1,2)){
         for(rep in 1:3){
           #Run the function
           factor.A <- 1
           factor.B <- 1
           factor.C <- 1
           factor.D <- 1 
           this.rep <- 1
           #Read in the R data file for a given setting
           load(paste(factor.A,".FactorA..",factor.B, ".FactorB..",
                 factor.C,".FactorC..",factor.D,".FactorD..",
                 this.rep, ".Rep.Rdata", sep = ""))
           setwd(home.dir)
           #Keep track of what setting you are on
           this.setting <-  paste(factor.A,".FactorA..",factor.B, ".FactorB..",
                                                 factor.C,".FactorC..",factor.D,".FactorD..",
                                                 this.rep, ".Rep", sep = "")

           #Run the master script; source it in
           source("Master_Analysis_Script_20240925.R")
           
           
           ##############################################################
           #### Append GWAS Results
           #Append results to the lists you created before the start of this rep
           these.spearman.rank.correlation.between.GWAS.core.QTNs[[count]] <-
             spearman.correlations.between.core.QTNs
           names(these.spearman.rank.correlation.between.GWAS.core.QTNs)[count] <-
               this.setting
           
           these.median.spearman.correlations.between.core.QTN[[count]] <-
             median.spearman.correlations.between.core.QTN
           names(these.spearman.rank.correlation.between.GWAS.core.QTNs)[count] <-
             this.setting
           
           
           these.median.core.spearman.correlation.within.types.of.selection.QTN[[count]] <-
             median.core.spearman.correlation.within.types.of.selection.QTN 
           names(these.median.core.spearman.correlation.within.types.of.selection.QTN)[count] <-
             this.setting
           
           these.median.core.spearman.correlation.different.types.of.selection.QTN[[count]] <-
             median.core.spearman.correlation.different.types.of.selection.QTN 
           names(these.median.core.spearman.correlation.different.types.of.selection.QTN)[count] <-
             this.setting
           
           
           
           these.spearman.rank.correlation.between.GWAS.peripheral.QTNs[[count]] <-
             spearman.correlations.between.peripheral.QTNs
           names(these.spearman.rank.correlation.between.GWAS.peripheral.QTNs)[count] <-
             this.setting 
           
           these.median.spearman.correlations.between.peripheral.QTN[[count]] <-
             median.spearman.correlations.between.peripheral.QTN
           names(these.median.spearman.correlations.between.peripheral.QTN)[count] <-
             this.setting
           
           these.median.peripheral.spearman.correlation.within.types.of.selection.QTN[[count]] <-
             median.peripheral.spearman.correlation.within.types.of.selection.QTN
           names(these.median.peripheral.spearman.correlation.within.types.of.selection.QTN)[count] <-
             this.setting
           
           
           these.median.peripheral.spearman.correlation.different.types.of.selection.QTN[[count]] <-
             median.peripheral.spearman.correlation.different.types.of.selection.QTN
           names(these.median.peripheral.spearman.correlation.different.types.of.selection.QTN)[count] <-
             this.setting
           
           these.spearman.rank.correlation.between.GWAS.core.SNPs[[count]] <-
             spearman.correlations.between.core.SNPs
           names(these.spearman.rank.correlation.between.GWAS.core.SNPs)[count] <-
             this.setting
           
           
           these.median.spearman.correlations.between.core.SNPs[[count]] <-
             median.spearman.correlations.between.core.SNPs
           names(these.median.spearman.correlations.between.core.SNPs)[count] <-
             this.setting
           
           
           these.median.core.spearman.correlation.within.types.of.selection.SNPs[[count]] <-
             median.core.spearman.correlation.within.types.of.selection.SNPs
           names(these.median.core.spearman.correlation.within.types.of.selection.SNPs)[count] <-
             this.setting
           
           these.median.core.spearman.correlation.different.types.of.selection.SNPs[[count]] <-
             median.core.spearman.correlation.different.types.of.selection.SNPs
           names(these.median.core.spearman.correlation.different.types.of.selection.SNPs)[count] <-
             this.setting
           
           
           these.spearman.rank.correlation.between.GWAS.peripheral.SNPs[[count]] <-
             spearman.correlations.between.peripheral.SNPs
           names(these.spearman.rank.correlation.between.GWAS.peripheral.SNPs)[count] <-
             this.setting  
           
           
           these.median.spearman.correlations.between.peripheral.SNPs[[count]] <-
             median.spearman.correlations.between.peripheral.SNPs
           names(these.median.spearman.correlations.between.peripheral.SNPs)[count] <-
             this.setting
           
           these.median.peripheral.spearman.correlation.within.types.of.selection.SNPs[[count]] <-
             median.peripheral.spearman.correlation.within.types.of.selection.SNPs
           names(these.median.peripheral.spearman.correlation.within.types.of.selection.SNPs)[count] <-
             this.setting
           
           these.median.peripheral.spearman.correlation.different.types.of.selection.SNPs[[count]] <-
             median.peripheral.spearman.correlation.different.types.of.selection.SNPs
           names(these.median.peripheral.spearman.correlation.different.types.of.selection.SNPs)[count] <-
             this.setting
           
   
           ##############################################################
           #### Append GS Results
           #Append results to the lists you created before the start of this rep
           these.prediction.accuracies.QTNs[[count]] <- prediction.accuracies.QTNs 
           names(these.prediction.accuracies.QTNs)[count] <- this.setting 
           
           these.prediction.accuracies.SNPs[[count]] <- prediction.accuracies.SNPs
           names(these.prediction.accuracies.SNPs)[count] <- this.setting 
           #Update count for the next iteration of this loop
           count <- count + 1
         }#End for(rep in 1:3)
      }#End for(el in 1:4)
    }#End for(k in 1:4)
  }#End for(j in 1:4)
}#End for(i in 1:3)

sessionInfo() 

