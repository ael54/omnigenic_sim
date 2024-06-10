


#Set your working directory
setwd("/Users/alipka/Library/CloudStorage/Box-Box/Sabbatical_Roslin_Institute/R_workspace/Sabbatical_Project")
home.dir <- getwd()

#####Read in all of the packages that are necessary
#Read in prerequiste libaries for GAPIT
library('MASS')
library(gplots)
library(sommer)



#Load the R data object containing all results
load("/Users/alipka/Library/CloudStorage/Box-Box/IR-281/analysis_results/analysis.1.FactorA..1.FactorB..0.05.FactorC..0.5.FactorD(7).Rdata")

#######################################################
### Results for the core QTL
#######################################################
#Make a plot summarizing the GWAS results of core QTL
#Initiate the plot
 #For loop through the pairs of three subpopulations
  #For loop through the main effects of Factors A-D
    #For loop through the different levels of Factor i
        ####
        #Box plot of Spearman rank correlation coefficients of core QTL (Y-axis)
        ### Against factor levels (X-axis)
    #End for loop through the different levels of Factor i
  #End for loop through the main effect pairs of populations
#End the plot
 

#########
#Make a plot of the variance component estimates of trait 2
#Initiate the plot
  #For loop through the pairs of phenotypic variance; genetic variance; breeding value variance
    #For loop through the main effects of Factors A-D
      #For loop through the different levels of Factor i
       ####
        # Plot variance component for trait 2 (Y-axis)
       ### Against factor levels (X-axis)
       ### Color-coded by generation (col = )
      #End for loop through the different levels of Factor i
    #End for loop through the main effects of Factors A-D
  #End for loop through the pairs of phenotypic variance; genetic variance; breeding value variance
#End the plot



#######################################################
### Results for the peripheral QTL
#######################################################
#Make a plot summarizing the GWAS results of core QTL
#Initiate the plot
  #For loop through the pairs of three subpopulations
    #For loop through the main effects of Factors A-D
      #For loop through the different levels of Factor i
        ####
        #Box plot of Spearman rank correlation coefficients of peripheral QTL (Y-axis)
        ### Against factor levels (X-axis)
    #End for loop through the different levels of Factor i
  #End for loop through the main effect pairs of populations
#End the plot


#########
#Make a plot of the variance component estimates of trait 3
#Initiate the plot
  #For loop through the pairs of phenotypic variance; genetic variance; breeding value variance
    #For loop through the main effects of Factors A-D
      #For loop through the different levels of Factor i
        ####
        # Plot variance component for trait 3 (Y-axis)
        ### Against factor levels (X-axis)
        ### Color-coded by generation (col = )
      #For loop through the different levels of Factor i
    #End for loop through the main effects of Factors A-D
  #End for loop through the pairs of phenotypic variance; genetic variance; breeding value variance
#End the plot



#######################################################################
### Results for the interaction between core and peripheral QTL
#######################################################################

#########
#Make a plot of the variance component estimates of trait 4
#Initiate the plot
  #For loop through the pairs of phenotypic variance; genetic variance; breeding value variance
    #For loop through the main effects of Factors A-D
      #For loop through the different levels of Factor i
        ####
        # Plot variance component for trait 4 (Y-axis)
        ### Against factor levels (X-axis)
        ### Color-coded by generation (col = )
      #For loop through the different levels of Factor i
    #End for loop through the main effects of Factors A-D
  #End for loop through the pairs of phenotypic variance; genetic variance; breeding value variance
#End the plot


#######################################################################
### Results for the genomic prediction analysis
#######################################################################
#Make a plot of the improvement of prediciton accuarices after accounting for the omingenic model
#Initiate the plot
  #For loop through the model that considers epsitasis versus the one that does not
    #For loop through the main effects of Factors A-D
      #For loop through the different levels of Factor i
         ####
         # Plot the improvement in prediction accuracy over a standard GBLUP model
         ### Against factor levels (X-axis)
      #For loop through the different levels of Factor i
    #End for loop through the main effects of Factors A-D
  #End for loop through the model that considers epsitasis versus the one that does not
#End the plot




########################################################################
#Old code, which might be useful.
# Delete all code below once a script is written
########################################################################
#This count variable will help with appending to these lists storing
# GWAS and GS summary statistics
count <- 1

#Initialize objects that will store GWAS results
these.spearman.rank.correlation.between.GWAS.core.QTNs <- list()
these.spearman.rank.correlation.between.GWAS.peripheral.QTNs <- list()

these.spearman.rank.correlation.between.GWAS.core.SNPs <- list()
these.spearman.rank.correlation.between.GWAS.peripheral.SNPs <- list()

#Initiallize objects that will store GS results
these.prediction.accuracies.QTNs <- list()
these.prediction.accuracies.SNPs <- list()

###For loop through all settings of the factorial experiment
### The end products populated versions of the updated lists
#### where each item of each list corresponds to results
#### from a particular setting of the factorial experiment
### Not to mention output directories of simulated traits

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
           source("Master_Analysis_Script_20240315.R")
           
           
           ##############################################################
           #### Append GWAS Results
           #Append results to the lists you created before the start of this rep
           these.spearman.rank.correlation.between.GWAS.core.QTNs[[count]] <-
             spearman.correlations.between.core.QTNs
           names(these.spearman.rank.correlation.between.GWAS.core.QTNs)[count] <-
               this.setting

           these.spearman.rank.correlation.between.GWAS.peripheral.QTNs[[count]] <-
             spearman.correlations.between.peripheral.QTNs
           names(these.spearman.rank.correlation.between.GWAS.peripheral.QTNs)[count] <-
             this.setting           
           
           these.spearman.rank.correlation.between.GWAS.core.SNPs[[count]] <-
             spearman.correlations.between.core.SNPs
           names(these.spearman.rank.correlation.between.GWAS.core.SNPs)[count] <-
             this.setting
           
           these.spearman.rank.correlation.between.GWAS.peripheral.SNPs[[count]] <-
             spearman.correlations.between.peripheral.SNPs
           names(these.spearman.rank.correlation.between.GWAS.peripheral.SNPs)[count] <-
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

