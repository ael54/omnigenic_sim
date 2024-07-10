


#Set your working directory
setwd("/Users/alipka/Library/CloudStorage/Box-Box/Sabbatical_Roslin_Institute/R_workspace/Sabbatical_Project")
home.dir <- getwd()

#####Read in all of the packages that are necessary
#Read in prerequiste libaries for GAPIT
library('MASS')
library(gplots)
library(sommer)
library(reshape)
library(ggplot2)


#######################################################
### Results for the core QTL
#######################################################
these.spearman.rank.correlation.between.GWAS.core.QTNs <- readRDS("/Users/alipka/Library/CloudStorage/Box-Box/IR-281/Master-Result-Lists-2024-06-13/master.these.spearman.rank.correlation.between.GWAS.core.QTNs.RDS")

#Extract all of the results and put them into an object
spearman.rank.vector <- NULL
factor.A.vector <- NULL
factor.B.vector <- NULL
factor.C.vector <- NULL
factor.D.vector <- NULL
rep.vector <- NULL
pop.A.vector <- NULL
pop.B.vector <- NULL
for(i in c(1,2,4)){
  for(j in c(0.05,0.5,1,2)){
    for(k in c(0.05,0.5,1,2)){
      for(el in c(0.05,0.5,1,2)){
        for(rep in 1:3){
          #Get the object from the list you want
          factor.A <- i
          factor.B <- j
          factor.C <- k
          factor.D <- el 
          this.rep <- rep
          
          #Read in the R data file for a given setting
          this.setting <-  paste(factor.A,".FactorA..",factor.B, ".FactorB..",
                                 factor.C,".FactorC..",factor.D,".FactorD..",
                                 this.rep, ".Rep", sep = "")
          
          #Extract the results we want
          results.this.setting <- these.spearman.rank.correlation.between.GWAS.core.QTNs[which(names(these.spearman.rank.correlation.between.GWAS.core.QTNs)
                                                                                               == this.setting)]
          
          
          #Put the results into a format that we can extract the numbers from                                                                                                                                                                        == this.setting)]
          results.this.setting.for.figure <- matrix(unlist(results.this.setting), nrow = 3)
          #Loop through the pairs of populations, and extract the information we need
          for(this.row in 1:nrow(results.this.setting.for.figure)){
            spearman.rank.vector <- c(spearman.rank.vector, results.this.setting.for.figure[this.row,3])
            factor.A.vector <- c(factor.A.vector, factor.A)
            factor.B.vector <- c(factor.B.vector, factor.B)
            factor.C.vector <- c(factor.C.vector, factor.C)
            factor.D.vector <- c(factor.D.vector, factor.D)
            rep.vector <- c(rep.vector, this.rep)
            pop.A.vector <-c(pop.A.vector, results.this.setting.for.figure[this.row,1])
            pop.B.vector <- c(pop.B.vector, results.this.setting.for.figure[this.row,2])
          }#End for(this.row in 1:nrow(results.this.setting.for.figure))
        }#End for(rep in 1:3)
      }#End for(el in 1:4)
    }#End for(k in 1:4)
  }#End for(j in 1:4)
}#End for(i in 1:3)

#Obtain the data for your boxplot
data.for.boxplot <- data.frame(spearman.rank.vector, factor.A.vector,
                               factor.B.vector, factor.C.vector,
                               factor.D.vector, rep.vector,
                               pop.A.vector, pop.B.vector)

#Make a plot summarizing the GWAS results of core QTL
#Ultimately, the information that I need for the box plots are
# Value: in column 1 and Factor Level in column 2

#Initiate the plot
pdf("GWAS.Results.Core.QTL.written.while.writing.plot.code.pdf", width = 8)
par(mfrow = c(3,4))
 #Two for loops through the pairs of three subpopulations
  for(this.pop.a in 1:2){
    for(this.pop.b in (this.pop.a+1):3){
        #Extract the rows for only the populations you want to compare
        data.for.boxplot.these.two.pops <- data.for.boxplot[which((data.for.boxplot$pop.A.vector == this.pop.a)
                                                                  &(data.for.boxplot$pop.B.vector == this.pop.b) ),]
        #Remove "NAs" from the results
        data.for.boxplot.these.two.pops <- na.omit(data.for.boxplot.these.two.pops)
        
      
        
        #Make a box plot for Factor A   
        means <- by(data.for.boxplot.these.two.pops$spearman.rank.vector, 
                    data.for.boxplot.these.two.pops$factor.A.vector, mean)
        
        boxplot(data.for.boxplot.these.two.pops$spearman.rank.vector ~ 
                data.for.boxplot.these.two.pops$factor.A.vector, col = "Red", ylim = c(-1,1),
                xlab = "Cor Add vs Per Add", ylab = paste("Spearman rank between populations ",
                                                         this.pop.a, " and ", this.pop.b, sep = ""))
        points(c(1:length(means)), means, pch = 3, cex = 0.75)
  
        
        #Make a box plot for Factor B
        means <- by(data.for.boxplot.these.two.pops$spearman.rank.vector, 
                    data.for.boxplot.these.two.pops$factor.B.vector, mean)
        
        boxplot(data.for.boxplot.these.two.pops$spearman.rank.vector ~ 
                  data.for.boxplot.these.two.pops$factor.B.vector, col = "Red", ylim = c(-1,1), 
                xlab = "Cor Epi vs Cor Add", ylab = paste("Spearman rank between populations ",
                                                         this.pop.a, " and ", this.pop.b, sep = ""))
        points(c(1:length(means)), means, pch = 3, cex = 0.75)
  
        
        #Make a box plot for Factor C
        means <- by(data.for.boxplot.these.two.pops$spearman.rank.vector, 
                    data.for.boxplot.these.two.pops$factor.C.vector, mean)
        
        boxplot(data.for.boxplot.these.two.pops$spearman.rank.vector ~ 
                  data.for.boxplot.these.two.pops$factor.C.vector, col = "Red", ylim = c(-1,1), 
                xlab = "Per Epi vs Per Add", ylab = paste("Spearman rank between populations ",
                                                         this.pop.a, " and ", this.pop.b, sep = ""))
        points(c(1:length(means)), means, pch = 3, cex = 0.75)
        
        
        #Make a box plot for Factor D
        means <- by(data.for.boxplot.these.two.pops$spearman.rank.vector, 
                    data.for.boxplot.these.two.pops$factor.D.vector, mean)
        
        boxplot(data.for.boxplot.these.two.pops$spearman.rank.vector ~ 
                  data.for.boxplot.these.two.pops$factor.D.vector, col = "Red", ylim = c(-1,1), 
                xlab = "Btw Epi vs Per Add", ylab = paste("Spearman rank between populations ",
                                                         this.pop.a, " and ", this.pop.b, sep = ""))
        points(c(1:length(means)), means, pch = 3, cex = 0.75)
        

          #For loop through the different levels of Factor i
              ####
              #Box plot of Spearman rank correlation coefficients of core QTL (Y-axis)
              ### Against factor levels (X-axis)
          #End for loop through the different levels of Factor i
      
        # End for loop through the main effect pairs of populations 
    }#End for(this.pop.b in (this.pop.a+1):3)
  }#End for(this.pop.a in 1:2) 
#End the plot
dev.off()

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
these.spearman.rank.correlation.between.GWAS.peripheral.QTNs  <- readRDS("/Users/alipka/Library/CloudStorage/Box-Box/IR-281/Master-Result-Lists-2024-06-13/master.these.spearman.rank.correlation.between.GWAS.peripheral.QTNs.RDS")

#Extract all of the results and put them into an object
spearman.rank.vector <- NULL
factor.A.vector <- NULL
factor.B.vector <- NULL
factor.C.vector <- NULL
factor.D.vector <- NULL
rep.vector <- NULL
pop.A.vector <- NULL
pop.B.vector <- NULL
for(i in c(1,2,4)){
  for(j in c(0.05,0.5,1,2)){
    for(k in c(0.05,0.5,1,2)){
      for(el in c(0.05,0.5,1,2)){
        for(rep in 1:3){
          #Get the object from the list you want
          factor.A <- i
          factor.B <- j
          factor.C <- k
          factor.D <- el 
          this.rep <- rep
          
          #Read in the R data file for a given setting
          this.setting <-  paste(factor.A,".FactorA..",factor.B, ".FactorB..",
                                 factor.C,".FactorC..",factor.D,".FactorD..",
                                 this.rep, ".Rep", sep = "")
          
          #Extract the results we want
          results.this.setting <- these.spearman.rank.correlation.between.GWAS.peripheral.QTNs[which(names(these.spearman.rank.correlation.between.GWAS.peripheral.QTNs)
                                                                                               == this.setting)]
          
          
          #Put the results into a format that we can extract the numbers from                                                                                                                                                                        == this.setting)]
          results.this.setting.for.figure <- matrix(unlist(results.this.setting), nrow = 3)
          #Loop through the pairs of populations, and extract the information we need
          for(this.row in 1:nrow(results.this.setting.for.figure)){
            spearman.rank.vector <- c(spearman.rank.vector, results.this.setting.for.figure[this.row,3])
            factor.A.vector <- c(factor.A.vector, factor.A)
            factor.B.vector <- c(factor.B.vector, factor.B)
            factor.C.vector <- c(factor.C.vector, factor.C)
            factor.D.vector <- c(factor.D.vector, factor.D)
            rep.vector <- c(rep.vector, this.rep)
            pop.A.vector <-c(pop.A.vector, results.this.setting.for.figure[this.row,1])
            pop.B.vector <- c(pop.B.vector, results.this.setting.for.figure[this.row,2])
          }#End for(this.row in 1:nrow(results.this.setting.for.figure))
        }#End for(rep in 1:3)
      }#End for(el in 1:4)
    }#End for(k in 1:4)
  }#End for(j in 1:4)
}#End for(i in 1:3)

#Obtain the data for your boxplot
data.for.boxplot <- data.frame(spearman.rank.vector, factor.A.vector,
                               factor.B.vector, factor.C.vector,
                               factor.D.vector, rep.vector,
                               pop.A.vector, pop.B.vector)

#Make a plot summarizing the GWAS results of core QTL
#Ultimately, the information that I need for the box plots are
# Value: in column 1 and Factor Level in column 2

#Initiate the plot
pdf("GWAS.Results.Peripheral.QTL.written.while.writing.plot.code.pdf", width = 8)
par(mfrow = c(3,4))
#Two for loops through the pairs of three subpopulations
for(this.pop.a in 1:2){
  for(this.pop.b in (this.pop.a+1):3){
    #Extract the rows for only the populations you want to compare
    data.for.boxplot.these.two.pops <- data.for.boxplot[which((data.for.boxplot$pop.A.vector == this.pop.a)
                                                              &(data.for.boxplot$pop.B.vector == this.pop.b) ),]
    #Remove "NAs" from the results
    data.for.boxplot.these.two.pops <- na.omit(data.for.boxplot.these.two.pops)
    
    
    
    #Make a box plot for Factor A   
    means <- by(data.for.boxplot.these.two.pops$spearman.rank.vector, 
                data.for.boxplot.these.two.pops$factor.A.vector, mean)
    
    boxplot(data.for.boxplot.these.two.pops$spearman.rank.vector ~ 
              data.for.boxplot.these.two.pops$factor.A.vector, col = "Blue", ylim = c(-1,1),
            xlab = "Cor Add vs Per Add", ylab = paste("Spearman rank between populations ",
                                                     this.pop.a, " and ", this.pop.b, sep = ""))
    points(c(1:length(means)), means, pch = 3, cex = 0.75)
    
    
    #Make a box plot for Factor B
    means <- by(data.for.boxplot.these.two.pops$spearman.rank.vector, 
                data.for.boxplot.these.two.pops$factor.B.vector, mean)
    
    boxplot(data.for.boxplot.these.two.pops$spearman.rank.vector ~ 
              data.for.boxplot.these.two.pops$factor.B.vector, col = "Blue", ylim = c(-1,1),
            xlab = "Cor Epi vs Cor Add", ylab = paste("Spearman rank between populations ",
                                                     this.pop.a, " and ", this.pop.b, sep = ""))
    points(c(1:length(means)), means, pch = 3, cex = 0.75)
    
    
    #Make a box plot for Factor C
    means <- by(data.for.boxplot.these.two.pops$spearman.rank.vector, 
                data.for.boxplot.these.two.pops$factor.C.vector, mean)
    
    boxplot(data.for.boxplot.these.two.pops$spearman.rank.vector ~ 
              data.for.boxplot.these.two.pops$factor.C.vector, col = "Blue", ylim = c(-1,1),
            xlab = "Per Epi vs Per Add", ylab = paste("Spearman rank between populations ",
                                                     this.pop.a, " and ", this.pop.b, sep = ""))
    points(c(1:length(means)), means, pch = 3, cex = 0.75)
    
    
    #Make a box plot for Factor D
    means <- by(data.for.boxplot.these.two.pops$spearman.rank.vector, 
                data.for.boxplot.these.two.pops$factor.D.vector, mean)
    
    boxplot(data.for.boxplot.these.two.pops$spearman.rank.vector ~ 
              data.for.boxplot.these.two.pops$factor.D.vector, col = "Blue", ylim = c(-1,1),
            xlab = "Btw Epi vs Per Add", ylab = paste("Spearman rank between populations ",
                                                     this.pop.a, " and ", this.pop.b, sep = ""))
    points(c(1:length(means)), means, pch = 3, cex = 0.75)
    
    
    #For loop through the different levels of Factor i
    ####
    #Box plot of Spearman rank correlation coefficients of core QTL (Y-axis)
    ### Against factor levels (X-axis)
    #End for loop through the different levels of Factor i
    
    # End for loop through the main effect pairs of populations 
  }#End for(this.pop.b in (this.pop.a+1):3)
}#End for(this.pop.a in 1:2) 
#End the plot
dev.off()


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
#######################################################
#######################################################
#######################################################
#Experimental code - delete once a pipeline has been developed
# maybe consider putting this into a function
these.prediction.accuracies.QTNs <- readRDS("/Users/alipka/Library/CloudStorage/Box-Box/IR-281/Master-Result-Lists-2024-06-13/master.these.prediction.accuracies.QTNs.RDS")

#####Get all of the results into an object that can be used to 
# make the box plots you want
validation.set.vector <- NULL
training.set.1.boolean <- NULL
training.set.2.boolean <- NULL
predictive.ability.GBLUP <- NULL
predictive.ability.Multi.Kern.add <- NULL
predictive.ability.Multi.Kern.epi <- NULL
factor.A.vector <- NULL
factor.B.vector <- NULL
factor.C.vector <- NULL
factor.D.vector <- NULL
rep.vector <- NULL

for(i in c(1,2,4)){
  for(j in c(0.05,0.5,1,2)){
    for(k in c(0.05,0.5,1,2)){
      for(el in c(0.05,0.5,1,2)){
        for(rep in 1:3){
          #Get the object from the list you want
          factor.A <- i
          factor.B <- j
          factor.C <- k
          factor.D <- el 
          this.rep <- rep
          
          this.setting <-  paste(factor.A,".FactorA..",factor.B, ".FactorB..",
                                 factor.C,".FactorC..",factor.D,".FactorD..",
                                 this.rep, ".Rep", sep = "")
          
          #Extract the results we want
          results.this.setting <- these.prediction.accuracies.QTNs[which(names(these.prediction.accuracies.QTNs)
                                                                         == this.setting)]
          
          
          #Put the results into a format that we can extract the numbers from                                                                                                                                                                        == this.setting)]
          results.this.setting.for.figure <- matrix(unlist(results.this.setting),nrow = 9)
          
          for(pred.abil in 1:nrow(results.this.setting.for.figure)){
            validation.set.vector <- c(validation.set.vector, results.this.setting.for.figure[pred.abil, 1])
            training.set.1.boolean <- c(training.set.1.boolean, results.this.setting.for.figure[pred.abil, 2])
            training.set.2.boolean <- c(training.set.2.boolean, results.this.setting.for.figure[pred.abil, 3])
            predictive.ability.GBLUP <- c(predictive.ability.GBLUP, results.this.setting.for.figure[pred.abil, 4])
            predictive.ability.Multi.Kern.add <- c(predictive.ability.Multi.Kern.add, results.this.setting.for.figure[pred.abil, 5])
            predictive.ability.Multi.Kern.epi <- c(predictive.ability.Multi.Kern.epi, results.this.setting.for.figure[pred.abil, 6])
          }#End for(pred.abil in 1:nrow(results.this.setting.for.figure))
          
          factor.A.vector <- c(factor.A.vector, rep(factor.A, nrow(results.this.setting.for.figure)))
          factor.B.vector <- c(factor.B.vector, rep(factor.B, nrow(results.this.setting.for.figure)))
          factor.C.vector <- c(factor.C.vector, rep(factor.C, nrow(results.this.setting.for.figure)))
          factor.D.vector <- c(factor.D.vector, rep(factor.D, nrow(results.this.setting.for.figure)))
          rep.vector <- c(rep.vector, rep(this.rep, nrow(results.this.setting.for.figure)))
          
        }#End for(rep in 1:3)
      }#End for(el in 1:4)
    }#End for(k in 1:4)
  }#End for(j in 1:4)
}#End for(i in 1:3)

data.for.boxplot <- data.frame(validation.set.vector, training.set.1.boolean, 
                               training.set.2.boolean, predictive.ability.GBLUP,
                               predictive.ability.Multi.Kern.add, predictive.ability.Multi.Kern.epi,
                               factor.A.vector, factor.B.vector, factor.C.vector, factor.D.vector,
                               rep.vector)



#######################################################
#######################################################
#######################################################
#######################################################

######Make a lattice plot of prediction accuracies. Rows = different training populations
# populations. Columns = Factors. X-axis = Factor levels, Y-axis = Prediction accuracies
# Separate box plots (on each plot) for each model investigated
#Initiate the plot
#for loop through each page of plots (validation population)
all.three.pops <- 1:length(unique(data.for.boxplot$validation.set.vector))
pdf("Box.plot.of.Prediciton.Accuracies.Experimental.20240709.pdf", width = 50)
  for(val.pop in all.three.pops){
   #For loop through each row (training sets)
    #Source in code that will make each row of the plot
      #Input parameters = data set for boxplot, variable indicate what the validation popluation is,
                        #variable indicating which training set it is
    #Source in the code that will create separate data sets for each combination of
    # training sets.
     
    these.two.training.sets <- all.three.pops[which(all.three.pops != val.pop)] 
    title.label <- paste("Val pop = ", val.pop, sep = "")
    source("/Users/alipka/Library/CloudStorage/Box-Box/Sabbatical_Roslin_Institute/R_workspace/Sabbatical_Project/Functions_to_Make_Life_Easier/Format_the_data_for_Box_Plots.R")
    
    
    #########Make the first row of the plot
    #Source in the code below that will help make these plots
    input.box.plot.data <- data.for.boxplot.this.validation.pop.melted.first.pop.train
    y.axis.label <- paste("Pop ",these.two.training.sets[1], " TS", sep = "")
    this.min.y.axis <- min(data.for.boxplot.this.validation.pop.melted.sep.for.each.train.pop$value)
    this.max.y.axis <- max(data.for.boxplot.this.validation.pop.melted.sep.for.each.train.pop$value)
    #Source in some gplot code that will make the plots, and save them as objects
    source("/Users/alipka/Library/CloudStorage/Box-Box/Sabbatical_Roslin_Institute/R_workspace/Sabbatical_Project/Functions_to_Make_Life_Easier/Make_GS_Box_Plots_20240619.R")
    fa.row.1 <- fa
    fb.row.1 <- fb
    fc.row.1 <- fc
    fd.row.1 <- fd
    
    #########Make the second row of the plot
    #Source in the code below that will help make these plots
    input.box.plot.data <- data.for.boxplot.this.validation.pop.melted.second.pop.train
    y.axis.label <- paste("Pop ",these.two.training.sets[2], " TS", sep = "")
    title.label <- NA
    
    #Source in some gplot code that will make the plots, and save them as objects
    source("/Users/alipka/Library/CloudStorage/Box-Box/Sabbatical_Roslin_Institute/R_workspace/Sabbatical_Project/Functions_to_Make_Life_Easier/Make_GS_Box_Plots_20240619.R")
    fa.row.2 <- fa
    fb.row.2 <- fb
    fc.row.2 <- fc
    fd.row.2 <- fd
    
    #########Make the third row of the plot
    #Source in the code below that will help make these plots
    input.box.plot.data <- data.for.boxplot.this.validation.pop.melted.both.pop.train
    y.axis.label <- paste("Pop ",these.two.training.sets[1], "+", these.two.training.sets[2], " TS", sep = "") 
    
    #Source in some gplot code that will make the plots, and save them as objects
    source("/Users/alipka/Library/CloudStorage/Box-Box/Sabbatical_Roslin_Institute/R_workspace/Sabbatical_Project/Functions_to_Make_Life_Easier/Make_GS_Box_Plots_20240619.R")
    fa.row.3 <- fa
    fb.row.3 <- fb
    fc.row.3 <- fc
    fd.row.3 <- fd
    
    #put everything into a box plot
    
   
     print(plot_grid(fa.row.1, fb.row.1, fc.row.1, fd.row.1,
              fa.row.2, fb.row.2, fc.row.2, fd.row.2,
              fa.row.3, fb.row.3, fc.row.3, fd.row.3,
              nrow = 3, ncol = 4))
  
  
    
    #########Make the second row of the plot
  }#End (val.pop in 1:length(unique(data.for.boxplot$validation.set.vector)))
    #End the plot
dev.off()


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

