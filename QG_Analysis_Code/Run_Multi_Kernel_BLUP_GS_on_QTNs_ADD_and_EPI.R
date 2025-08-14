#Run a GS for each subpopulation
list.of.subpopulation.traits <- list(directional.subpopulation.trait.10.pct, directional.subpopulation.trait.20.pct, 
                                     disruptive.subpopulation.trait.10.pct, disruptive.subpopulation.trait.20.pct,
                                     stabilizing.subpopulation.trait.10.pct, stabilizing.subpopulation.trait.20.pct,
                                     directional.subpopulation.trait.10.prev.gen, directional.subpopulation.trait.20.prev.gen,
                                     disruptive.subpopulation.trait.10.prev.gen,disruptive.subpopulation.trait.20.prev.gen,
                                     stabilizing.subpopulation.trait.10.prev.gen,stabilizing.subpopulation.trait.20.prev.gen)

list.of.subpopulation.QTN <- list(directional.subpopulation.QTNs.10.pct, directional.subpopulation.QTNs.20.pct,
                                  disruptive.subpopulation.QTNs.10.pct, disruptive.subpopulation.QTNs.20.pct,
                                  stabilizing.subpopulation.QTNs.10.pct, stabilizing.subpopulation.QTNs.20.pct,
                                  directional.subpopulation.QTNs.10.prev.gen, directional.subpopulation.QTNs.20.prev.gen,
                                  disruptive.subpopulation.QTNs.10.prev.gen,disruptive.subpopulation.QTNs.20.prev.gen,
                                  stabilizing.subpopulation.QTNs.10.prev.gen, stabilizing.subpopulation.QTNs.20.prev.gen)

names.of.subpopulations <- c("Directional.selection.10.pct","Directional.selection.20.pct",
                             "Disruptive.selection.10.pct", "Disruptive.selection.20.pct",
                             "Stabilizing.selection.10.pct","Stabilizing.selection.20.pct",
                             "Directional.selection.10.pct.prev.gen","Directional.selection.20.pct.prev.gen",
                             "Disruptive.selection.10.pct.prev.gen", "Disruptive.selection.20.pct.prev.gen",
                             "Stabilizing.selection.10.pct.prev.gen","Stabilizing.selection.20.pct.prev.gen")


#Obtain the prediction accuracies using Multi_Kernel_BLUP
validation.set <- NULL
training.set.0 <- NULL
training.set.1 <- NULL
training.set.2 <- NULL
training.set.3a <- NULL
training.set.3b <- NULL
prediction.accuracy.Multi.Kern.No.Epi <- NULL
prediction.accuracy.Multi.Kern.With.Epi <- NULL
for(eye in 1:(length(names.of.subpopulations)/2)){
  pct.10 <- FALSE
  pct.20 <- FALSE
  
  #Choose a particular subpopulation to be the validation population
  this.validation.set <- eye
  if(length(grep("10", names.of.subpopulations[eye])) != 0) pct.10 = TRUE
  if(length(grep("20", names.of.subpopulations[eye])) != 0) pct.20 = TRUE 
  
  
  this.training.set.0.names <- names.of.subpopulations[which(substr(names.of.subpopulations,start = 1, stop = 3) 
                                                             == substr(names.of.subpopulations[eye], start = 1, stop = 3))]
  this.training.set.0.names <- this.training.set.0.names[grep(".prev.gen", this.training.set.0.names)]
  if(pct.10) this.training.set.0.names <- this.training.set.0.names[grep("10", this.training.set.0.names)]
  if(pct.20) this.training.set.0.names <- this.training.set.0.names[grep("20", this.training.set.0.names)]
  this.training.set.0 <- which(names.of.subpopulations == this.training.set.0.names)
  
  
  this.training.set.1.names <- names.of.subpopulations[which(substr(names.of.subpopulations,start = 1, stop = 3) 
                                                             == substr(names.of.subpopulations[eye], start = 1, stop = 3))]
  this.training.set.1.names <- this.training.set.1.names[-which(this.training.set.1.names  %in% names.of.subpopulations[eye])]
  this.training.set.1.names <- this.training.set.1.names[-grep(".prev.gen", this.training.set.1.names)]
  this.training.set.1 <- which(names.of.subpopulations == this.training.set.1.names)
  
  
  this.training.set.2.names <- names.of.subpopulations[-which(substr(names.of.subpopulations,start = 1, stop = 3) 
                                                              == substr(names.of.subpopulations[eye], start = 1, stop = 3))]
  this.training.set.2.names <- this.training.set.2.names[-grep(".prev.gen", this.training.set.2.names)]
  this.training.set.2 <- which(names.of.subpopulations %in% this.training.set.2.names)
  
  
  this.training.set.3a.names <- names.of.subpopulations[-eye]
  this.training.set.3a.names <- this.training.set.3a.names[-grep(".prev.gen", this.training.set.3a.names)]
  this.training.set.3a <- which(names.of.subpopulations %in% this.training.set.3a.names)
  
  this.training.set.3b.names <- c(this.training.set.3a.names, this.training.set.0.names)
  this.training.set.3b <- which(names.of.subpopulations %in% this.training.set.3b.names)
  
  
  ##################################################################
  #Run genomic prediction using one of the other two populations as 
  # the training set
  this.training.set <- this.training.set.0
  validation.set <-c(validation.set, this.validation.set)
  training.set.0 <- c(training.set.0, 1)  
  training.set.1 <- c(training.set.1, 0)
  training.set.2 <- c(training.set.2, 0)
  training.set.3a <- c(training.set.3a, 0)
  training.set.3b <- c(training.set.3b, 0)
  
  #####Multi-kernel model where the additive and epistatic effects of core and peripheral
  # QTN are separate
  source("QG_Analysis_Code/Run_Multi_Kernel_BLUP_Add_and_Epi.R")
  prediction.accuracy.Multi.Kern.With.Epi <- c(prediction.accuracy.Multi.Kern.With.Epi, r.gy.add.epi.mult.kern)
  print("Yes! Run_Multi_Kernel_BLUP_Add_and_Epi.R run successfully on QTNs!")
  
  # the training set
  this.training.set <- this.training.set.1
  validation.set <-c(validation.set, this.validation.set)
  training.set.0 <- c(training.set.0, 0)  
  training.set.1 <- c(training.set.1, 1)
  training.set.2 <- c(training.set.2, 0)
  training.set.3a <- c(training.set.3a, 0)
  training.set.3b <- c(training.set.3b, 0)
  
  #####Multi-kernel model where the additive and epistatic effects of core and peripheral
  # QTN are separate
  source("QG_Analysis_Code/Run_Multi_Kernel_BLUP_Add_and_Epi.R")
  prediction.accuracy.Multi.Kern.With.Epi <- c(prediction.accuracy.Multi.Kern.With.Epi, r.gy.add.epi.mult.kern)
  print("Yes! Run_Multi_Kernel_BLUP_Add_and_Epi.R run successfully on QTNs AGAIN!")
  
  ##################################################################
  #Run genomic prediction using the second of the other two populations as 
  # the training set
  this.training.set <- this.training.set.2
  validation.set <-c(validation.set, this.validation.set)
  training.set.0 <- c(training.set.0, 0)  
  training.set.1 <- c(training.set.1, 0)
  training.set.2 <- c(training.set.2, 1)
  training.set.3a <- c(training.set.3a, 0)
  training.set.3b <- c(training.set.3b, 0)
  
  #####Multi-kernel model where the additive and epistatic effects of core and peripheral
  # QTN are separate
  source("QG_Analysis_Code/Run_Multi_Kernel_BLUP_Add_and_Epi.R")
  prediction.accuracy.Multi.Kern.With.Epi <- c(prediction.accuracy.Multi.Kern.With.Epi, r.gy.add.epi.mult.kern)
  print("Yes! Run_Multi_Kernel_BLUP_Add_and_Epi.R run successfully on QTNs AGAIN 2!")
  
  ###############################################
  # the training set
  this.training.set <- this.training.set.3a
  validation.set <-c(validation.set, this.validation.set)
  training.set.0 <- c(training.set.0, 0)  
  training.set.1 <- c(training.set.1, 0)
  training.set.2 <- c(training.set.2, 0)
  training.set.3a <- c(training.set.3a, 1)
  training.set.3b <- c(training.set.3b, 0)
  
  #####Multi-kernel model where the additive and epistatic effects of core and peripheral
  # QTN are separate
  source("QG_Analysis_Code/Run_Multi_Kernel_BLUP_Add_and_Epi.R")
  prediction.accuracy.Multi.Kern.With.Epi <- c(prediction.accuracy.Multi.Kern.With.Epi, r.gy.add.epi.mult.kern)
  print("Yes! Run_Multi_Kernel_BLUP_Add_and_Epi.R run successfully on QTNs AGAIN 3!")
  #Standard GBLUP model
  ###############################################
  # the training set
  this.training.set <- this.training.set.3b
  validation.set <-c(validation.set, this.validation.set)
  training.set.0 <- c(training.set.0, 0)  
  training.set.1 <- c(training.set.1, 0)
  training.set.2 <- c(training.set.2, 0)
  training.set.3a <- c(training.set.3a, 0)
  training.set.3b <- c(training.set.3b, 1)
  
  #####Multi-kernel model where the additive and epistatic effects of core and peripheral
  # QTN are separate
  source("QG_Analysis_Code/Run_Multi_Kernel_BLUP_Add_and_Epi.R")
  prediction.accuracy.Multi.Kern.With.Epi <- c(prediction.accuracy.Multi.Kern.With.Epi, r.gy.add.epi.mult.kern)
  print("Yes! Run_Multi_Kernel_BLUP_Add_and_Epi.R run successfully on QTNs AGAIN 4!")
  
} # End for(i in 1:length(names.of.subpopulations))
prediction.accuracies.QTNs  <- data.frame(validation.set, training.set.0, training.set.1, 
                                          training.set.2, training.set.3a, training.set.3b,
                                          prediction.accuracy.Multi.Kern.With.Epi)



