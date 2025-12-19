
for (eye in rest_of_eye) {
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
  source("Special_Multi_Kernel_ADD_EPI_QTNs/Run_Multi_Kernel_BLUP_Add_and_Epi.R")
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
  source("Special_Multi_Kernel_ADD_EPI_QTNs/Run_Multi_Kernel_BLUP_Add_and_Epi.R")
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
  source("Special_Multi_Kernel_ADD_EPI_QTNs/Run_Multi_Kernel_BLUP_Add_and_Epi.R")
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
  source("Special_Multi_Kernel_ADD_EPI_QTNs/Run_Multi_Kernel_BLUP_Add_and_Epi.R")
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
  source("Special_Multi_Kernel_ADD_EPI_QTNs/Run_Multi_Kernel_BLUP_Add_and_Epi.R")
  prediction.accuracy.Multi.Kern.With.Epi <- c(prediction.accuracy.Multi.Kern.With.Epi, r.gy.add.epi.mult.kern)
  print("Yes! Run_Multi_Kernel_BLUP_Add_and_Epi.R run successfully on QTNs AGAIN 4!")
  
}

prediction.accuracies.QTNs  <- data.frame(validation.set, training.set.0, training.set.1, 
                                          training.set.2, training.set.3a, training.set.3b,
                                          prediction.accuracy.Multi.Kern.With.Epi)