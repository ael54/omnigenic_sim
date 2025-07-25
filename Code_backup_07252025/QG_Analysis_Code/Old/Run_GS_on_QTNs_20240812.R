#Run a GS for each subpopulation
list.of.subpopulation.traits <- list(directional.subpopulation.trait.10.pct, directional.subpopulation.trait.20.pct, 
                                     disruptive.subpopulation.trait.10.pct, disruptive.subpopulation.trait.20.pct,
                                     stabilizing.subpopulation.trait.10.pct, stabilizing.subpopulation.trait.20.pct)
list.of.subpopulation.QTN <- list(directional.subpopulation.QTNs.10.pct, directional.subpopulation.QTNs.20.pct,
                                  disruptive.subpopulation.QTNs.10.pct, disruptive.subpopulation.QTNs.20.pct,
                                  stabilizing.subpopulation.QTNs.10.pct, stabilizing.subpopulation.QTNs.20.pct)
#list.of.subpopulation.SNPs <- list(directional.subpopulation.SNPs,
#                                   disruptive.subpopulation.SNPs,
#                                   stabilizing.subpopulation.SNPs)
names.of.subpopulations <- c("Directional.selection.10.pct","Directional.selection.20.pct",
                             "Disruptive.selection.10.pct", "Disruptive.selection.20.pct",
                             "Stabilizing.selection.10.pct","Stabilizing.selection.20.pct")


#Obtain the prediction accuracies using GBLUP
validation.set <- NULL
training.set.1 <- NULL
training.set.2 <- NULL
training.set.3 <- NULL
prediction.accuracy.GBLUP <- NULL
for(eye in 1:length(names.of.subpopulations)){
  #Choose a particular subpopulation to be the validation population
  this.validation.set <- eye

  this.training.set.1.names <- names.of.subpopulations[which(substr(names.of.subpopulations,start = 1, stop = 3) 
                               == substr(names.of.subpopulations[eye], start = 1, stop = 3))]
  this.training.set.1.names <- this.training.set.1.names[-which(this.training.set.1.names  %in% names.of.subpopulations[eye])]
  this.training.set.1 <- which(names.of.subpopulations == this.training.set.1.names)
  
  
  this.training.set.2.names <- names.of.subpopulations[-which(substr(names.of.subpopulations,start = 1, stop = 3) 
                                                       == substr(names.of.subpopulations[eye], start = 1, stop = 3))]
  this.training.set.2 <- which(names.of.subpopulations %in% this.training.set.2.names)
  
  
  this.training.set.3.names <- names.of.subpopulations[-eye]
  this.training.set.3 <- which(names.of.subpopulations %in% this.training.set.3.names)
  
  
  ##################################################################
  #Run genomic prediction using one of the other two populations as 
  # the training set
  this.training.set <- this.training.set.1
  validation.set <-c(validation.set, this.validation.set)
  training.set.1 <- c(training.set.1, 1)
  training.set.2 <- c(training.set.2, 0)
  training.set.3 <- c(training.set.3, 0)
  
  ####Standard GBLUP model
  source("QG_Analysis_Code/Run_GBLUP_Train_in_One_Pop_Validate_in_Another.R")
  prediction.accuracy.GBLUP <- c(prediction.accuracy.GBLUP, r.gy) 
 
 
  ##################################################################
  #Run genomic prediction using the second of the other two populations as 
  # the training set
  this.training.set <- this.training.set.2
  validation.set <-c(validation.set, this.validation.set)
  training.set.1 <- c(training.set.1, 0)
  training.set.2 <- c(training.set.2, 1)
  training.set.3 <- c(training.set.3, 0)
  
  #Standard GBLUP model
  source("QG_Analysis_Code/Run_GBLUP_Train_in_One_Pop_Validate_in_Another.R")
  prediction.accuracy.GBLUP <- c(prediction.accuracy.GBLUP, r.gy) 
  
  
  # the training set
  this.training.set <- this.training.set.3
  validation.set <-c(validation.set, this.validation.set)
  training.set.1 <- c(training.set.1, 0)
  training.set.2 <- c(training.set.2, 0)
  training.set.3 <- c(training.set.3, 1)
  
  #Standard GBLUP model
  source("QG_Analysis_Code/Run_GBLUP_Train_in_One_Pop_Validate_in_Another.R")
  prediction.accuracy.GBLUP <- c(prediction.accuracy.GBLUP, r.gy) 
  
  
  #Multi-kernel model where the additive effects of core and peripheral
  # QTN are separate
  
  #Multi-kernel model where the additive and epistatic effects of core and peripheral
  # QTN are separate
  
  #
  
  
  } # End for(i in 1:length(names.of.subpopulations))

#To-do: compile the results into a data frame
prediction.accuracies.QTNs <- data.frame(validation.set, training.set.1, training.set.2 ,
                                    prediction.accuracy.GBLUP)

#Obtain the prediction accuracies for the multi-kernel-models
validation.set <- NULL
training.set.1 <- NULL
training.set.2 <- NULL
prediction.accuracy.Multi.Kern.No.Epi <- NULL
prediction.accuracy.Multi.Kern.With.Epi <- NULL
for(eye in 1:length(names.of.subpopulations)){
  #Choose a particular subpopulation to be the validation population
  #Choose a particular subpopulation to be the validation population
  this.validation.set <- eye
  
  this.training.set.1.names <- names.of.subpopulations[which(substr(names.of.subpopulations,start = 1, stop = 3) 
                                                             == substr(names.of.subpopulations[eye], start = 1, stop = 3))]
  this.training.set.1.names <- this.training.set.1.names[-which(this.training.set.1.names  %in% names.of.subpopulations[eye])]
  this.training.set.1 <- which(names.of.subpopulations == this.training.set.1.names)
  
  
  this.training.set.2.names <- names.of.subpopulations[-which(substr(names.of.subpopulations,start = 1, stop = 3) 
                                                              == substr(names.of.subpopulations[eye], start = 1, stop = 3))]
  this.training.set.2 <- which(names.of.subpopulations %in% this.training.set.2.names)
  
  
  this.training.set.3.names <- names.of.subpopulations[-eye]
  this.training.set.3 <- which(names.of.subpopulations %in% this.training.set.3.names)
  
  
  
  ##################################################################
  #Run genomic prediction using one of the other two populations as 
  # the training set
  this.training.set <- this.training.set.1
  validation.set <-c(validation.set, this.validation.set)
  training.set.1 <- c(training.set.1, 1)
  training.set.2 <- c(training.set.2, 0)
  training.set.3 <- c(training.set.3, 0)
  
  #####Multi-kernel model where the additive effects of core and peripheral
  # QTN are separate
  source("QG_Analysis_Code/Run_Multi_Kernel_BLUP_Add_Only.R")
  prediction.accuracy.Multi.Kern.No.Epi <- c(prediction.accuracy.Multi.Kern.No.Epi,  r.gy.add.mult.kern)
  
  #####Multi-kernel model where the additive and epistatic effects of core and peripheral
  # QTN are separate
  source("QG_Analysis_Code/Run_Multi_Kernel_BLUP_Add_and_Epi.R")
  prediction.accuracy.Multi.Kern.With.Epi <- c(prediction.accuracy.Multi.Kern.With.Epi, r.gy.add.epi.mult.kern)
  ##################################################################
  #Run genomic prediction using the second of the other two populations as 
  # the training set
  this.training.set <- this.training.set.2
  validation.set <-c(validation.set, this.validation.set)
  training.set.1 <- c(training.set.1, 0)
  training.set.2 <- c(training.set.2, 1)
  training.set.3 <- c(training.set.3, 0)
  
  
  
  #####Multi-kernel model where the additive effects of core and peripheral
  # QTN are separate
  source("QG_Analysis_Code/Run_Multi_Kernel_BLUP_Add_Only.R")
  prediction.accuracy.Multi.Kern.No.Epi <- c(prediction.accuracy.Multi.Kern.No.Epi,  r.gy.add.mult.kern)
  
  #####Multi-kernel model where the additive and epistatic effects of core and peripheral
  # QTN are separate
  source("QG_Analysis_Code/Run_Multi_Kernel_BLUP_Add_and_Epi.R")
  prediction.accuracy.Multi.Kern.With.Epi <- c(prediction.accuracy.Multi.Kern.With.Epi, r.gy.add.epi.mult.kern)
  
  ##################################################################
  #Run genomic prediction using the second of both of the other two 
  # populations as the training set
  this.training.set <- this.training.set.3
  validation.set <-c(validation.set, this.validation.set)
  training.set.1 <- c(training.set.1, 0)
  training.set.2 <- c(training.set.2, 0)
  training.set.3 <- c(training.set.3, 1)
  
   
  
  #####Multi-kernel model where the additive effects of core and peripheral
  # QTN are separate
  source("QG_Analysis_Code/Run_Multi_Kernel_BLUP_Add_Only.R")
  prediction.accuracy.Multi.Kern.No.Epi <- c(prediction.accuracy.Multi.Kern.No.Epi,  r.gy.add.mult.kern)
  
  #####Multi-kernel model where the additive and epistatic effects of core and peripheral
  # QTN are separate
  source("QG_Analysis_Code/Run_Multi_Kernel_BLUP_Add_and_Epi.R")
  prediction.accuracy.Multi.Kern.With.Epi <- c(prediction.accuracy.Multi.Kern.With.Epi, r.gy.add.epi.mult.kern)
  
} # End for(i in 1:length(names.of.subpopulations))

prediction.accuracies.QTNs  <- data.frame(prediction.accuracies.QTNs, prediction.accuracy.Multi.Kern.No.Epi,
                                    prediction.accuracy.Multi.Kern.With.Epi)



