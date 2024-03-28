#Run a GWAS using a GLM for each subpopulation
list.of.subpopulation.traits <- list(directional.subpopulation.trait,
                                     disruptive.subpopulation.trait,
                                     stabilizing.subpopulation.trait)
#list.of.subpopulation.QTN <- list(directional.subpopulation.QTNs,
#                                  disruptive.subpopulation.QTNs,
#                                  stabilizing.subpopulation.QTNs)
list.of.subpopulation.SNPs <- list(directional.subpopulation.SNPs,
                                   disruptive.subpopulation.SNPs,
                                   stabilizing.subpopulation.SNPs)
names.of.subpopulations <- c("Directional.selection", "Disruptive.selection",
                             "Stabilizing.selection")


#Obtain the prediction accuracies using GBLUP
validation.set <- NULL
training.set.1 <- NULL
training.set.2 <- NULL
prediction.accuracy.GBLUP <- NULL
for(i in 1:length(names.of.subpopulations)){
  #Choose a particular subpopulation to be the validation population
  this.validation.set <- i
  this.training.set.1 <- which(1:length(names.of.subpopulations)!= i)[1]
  this.training.set.2 <- which(1:length(names.of.subpopulations)!= i)[2]
  
  ##################################################################
  #Run genomic prediction using one of the other two populations as 
  # the training set
  this.training.set <- this.training.set.1
  validation.set <-c(validation.set, this.validation.set)
  training.set.1 <- c(training.set.1, 1)
  training.set.2 <- c(training.set.2, 0)
  ####Standard GBLUP model
  source("QG_Analysis_Code/Run_GBLUP_Train_in_One_Pop_Validate_in_Another_using_SNPs.R")
  prediction.accuracy.GBLUP <- c(prediction.accuracy.GBLUP, r.gy) 
 
 
  ##################################################################
  #Run genomic prediction using the second of the other two populations as 
  # the training set
  this.training.set <- this.training.set.2
  validation.set <-c(validation.set, this.validation.set)
  training.set.1 <- c(training.set.1, 0)
  training.set.2 <- c(training.set.2, 1)
  
  #Standard GBLUP model
  source("QG_Analysis_Code/Run_GBLUP_Train_in_One_Pop_Validate_in_Another_using_SNPs.R")
  prediction.accuracy.GBLUP <- c(prediction.accuracy.GBLUP, r.gy) 
  
  
  #Multi-kernel model where the additive effects of core and peripheral
  # QTN are separate
  
  #Multi-kernel model where the additive and epistatic effects of core and peripheral
  # QTN are separate
  
  ##################################################################
  #Run genomic prediction using the second of both of the other two 
  # populations as the training set
  this.training.set <- c(this.training.set.1, this.training.set.2)
  validation.set <-c(validation.set, this.validation.set)
  training.set.1 <- c(training.set.1, 1)
  training.set.2 <- c(training.set.2, 1)
  
  #Standard GBLUP model
  source("QG_Analysis_Code/Run_GBLUP_Train_in_One_Pop_Validate_in_Another_using_SNPs.R")
  prediction.accuracy.GBLUP <- c(prediction.accuracy.GBLUP, r.gy) 
  
  
  #Multi-kernel model where the additive effects of core and peripheral
  # QTN are separate
  
  #Multi-kernel model where the additive and epistatic effects of core and peripheral
  # QTN are separate
  
  
  
  
  
  } # End for(i in 1:length(names.of.subpopulations))

#To-do: compile the results into a data frame
prediction.accuracies.SNPs <- data.frame(validation.set, training.set.1, training.set.2 ,
                                    prediction.accuracy.GBLUP)
library(sommer)
#Obtain the prediction accuracies for the multi-kernel-models
validation.set <- NULL
training.set.1 <- NULL
training.set.2 <- NULL
prediction.accuracy.Multi.Kern.No.Epi <- NULL
prediction.accuracy.Multi.Kern.With.Epi <- NULL
for(i in 1:length(names.of.subpopulations)){
  #Choose a particular subpopulation to be the validation population
  this.validation.set <- i
  this.training.set.1 <- which(1:length(names.of.subpopulations)!= i)[1]
  this.training.set.2 <- which(1:length(names.of.subpopulations)!= i)[2]
  
  ##################################################################
  #Run genomic prediction using one of the other two populations as 
  # the training set
  this.training.set <- this.training.set.1
  validation.set <-c(validation.set, this.validation.set)
  training.set.1 <- c(training.set.1, 1)
  training.set.2 <- c(training.set.2, 0)
  
  #####Multi-kernel model where the additive effects of core and peripheral
  # QTN are separate
  source("QG_Analysis_Code/Run_Multi_Kernel_BLUP_Add_Only_using_SNPs.R")
  prediction.accuracy.Multi.Kern.No.Epi <- c(prediction.accuracy.Multi.Kern.No.Epi,  r.gy.add.mult.kern)
  
  #####Multi-kernel model where the additive and epistatic effects of core and peripheral
  # QTN are separate
  source("QG_Analysis_Code/Run_Multi_Kernel_BLUP_Add_and_Epi_using_SNPs.R")
  prediction.accuracy.Multi.Kern.With.Epi <- c(prediction.accuracy.Multi.Kern.With.Epi, r.gy.add.epi.mult.kern)
  ##################################################################
  #Run genomic prediction using the second of the other two populations as 
  # the training set
  this.training.set <- this.training.set.2
  validation.set <-c(validation.set, this.validation.set)
  training.set.1 <- c(training.set.1, 0)
  training.set.2 <- c(training.set.2, 1)
  
  
  
  #####Multi-kernel model where the additive effects of core and peripheral
  # QTN are separate
  source("QG_Analysis_Code/Run_Multi_Kernel_BLUP_Add_Only_using_SNPs.R")
  prediction.accuracy.Multi.Kern.No.Epi <- c(prediction.accuracy.Multi.Kern.No.Epi,  r.gy.add.mult.kern)
  
  #####Multi-kernel model where the additive and epistatic effects of core and peripheral
  # QTN are separate
  source("QG_Analysis_Code/Run_Multi_Kernel_BLUP_Add_and_Epi_using_SNPs.R")
  prediction.accuracy.Multi.Kern.With.Epi <- c(prediction.accuracy.Multi.Kern.With.Epi, r.gy.add.epi.mult.kern)
  
  ##################################################################
  #Run genomic prediction using the second of both of the other two 
  # populations as the training set
  this.training.set <- c(this.training.set.1, this.training.set.2)
  validation.set <-c(validation.set, this.validation.set)
  training.set.1 <- c(training.set.1, 1)
  training.set.2 <- c(training.set.2, 1)
  
   
  
  #####Multi-kernel model where the additive effects of core and peripheral
  # QTN are separate
  source("QG_Analysis_Code/Run_Multi_Kernel_BLUP_Add_Only_using_SNPs.R")
  prediction.accuracy.Multi.Kern.No.Epi <- c(prediction.accuracy.Multi.Kern.No.Epi,  r.gy.add.mult.kern)
  
  #####Multi-kernel model where the additive and epistatic effects of core and peripheral
  # QTN are separate
  source("QG_Analysis_Code/Run_Multi_Kernel_BLUP_Add_and_Epi_using_SNPs.R")
  prediction.accuracy.Multi.Kern.With.Epi <- c(prediction.accuracy.Multi.Kern.With.Epi, r.gy.add.epi.mult.kern)
  
} # End for(i in 1:length(names.of.subpopulations))

prediction.accuracies.SNPs <- data.frame(prediction.accuracies.SNPs, prediction.accuracy.Multi.Kern.No.Epi,
                                    prediction.accuracy.Multi.Kern.With.Epi)



