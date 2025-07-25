


#Set your working directory
setwd("/Users/alipka/Library/CloudStorage/Box-Box/Sabbatical_Roslin_Institute/R_workspace/omnigenic_sim")
home.dir <- getwd()

#####Read in all of the packages that are necessary
#Read in prerequiste libaries for GAPIT
library('MASS')
library(gplots)
library(sommer)
library(reshape)
library(ggplot2)
library(cowplot)



#########
#Make a plot of the variance component estimates of trait 2
#Read in the results for the trait values, genetic values, and breeding values
vc.traits <- readRDS("/Users/alipka/Library/CloudStorage/Box-Box/IR-281/prelim-results-analysis-pipeline-20241030/varcovar-analysis-pipeline-20241030/master.these.breeding.value.var.covar.RDS") 
vc.genetic.values <- readRDS("/Users/alipka/Library/CloudStorage/Box-Box/IR-281/prelim-results-analysis-pipeline-20241030/varcovar-analysis-pipeline-20241030/master.these.genetic.value.var.covar.RDS")
vc.breeding.values <- readRDS("/Users/alipka/Library/CloudStorage/Box-Box/IR-281/prelim-results-analysis-pipeline-20241030/varcovar-analysis-pipeline-20241030/master.these.trait.var.covar.RDS")
#Format the results for trait values, genetic values, and breeding values so that they can be
# used to make graphs



#Create three objects: one for trait values, one for genetic values, and one for breeding values
# Each of these objects will have information for traits 1-4 for each generation
trait.1.variance <- NULL
trait.2.variance <- NULL
trait.3.variance <- NULL
trait.4.variance <- NULL
genetic.value.1.variance <- NULL
genetic.value.2.variance <- NULL
genetic.value.3.variance <- NULL
genetic.value.4.variance <- NULL
breeding.value.1.variance <- NULL
breeding.value.2.variance <- NULL
breeding.value.3.variance <- NULL
breeding.value.4.variance <- NULL
factor.A.vector <- NULL
factor.B.vector <- NULL
factor.C.vector <- NULL
factor.D.vector <- NULL
rep.vector <- NULL
subpopulation.vector <- NULL
selection.number.vector <- NULL
generation.vector <- NULL
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
            this.setting <-  paste(".Factor.A.", i,
                                   "Factor.B.", j,
                                   "Factor.C.", k,
                                   "Factor.D.", el,
                                   "Rep.", this.rep,
                                   sep = "")
            
            #Extract the results we want
            vc.traits.this.setting <- vc.traits[which(grepl(this.setting, names(vc.traits)))]
            vc.genetic.values.this.setting <- vc.genetic.values[which(grepl(this.setting, names(vc.genetic.values)))]  
            vc.breeding.values.this.setting <- vc.breeding.values[which(grepl(this.setting, names(vc.breeding.values)))]  
            
            #Move on to the next iteration if thre are no result
            if(length(vc.traits.this.setting) == 0){
              print(paste("Aw snap! There are no results for ", this.setting, sep = ""))
              next
            }# End if(length(vc.traits.this.setting) == 0)
            
            #Initiate everything for the founding population
            factor.A.vector <- c(factor.A.vector, factor.A)
            factor.B.vector <- c(factor.B.vector, factor.B)
            factor.C.vector <- c(factor.C.vector, factor.C)
            factor.D.vector <- c(factor.D.vector, factor.D)
            rep.vector <- c(rep.vector, this.rep)  
            subpopulation.vector <- c(subpopulation.vector, "Founder")
            selection.number.vector <- c(selection.number.vector, "Founder")
            generation.vector <- c(generation.vector, 0)
            
            
            this.trait.var.covar <- matrix(unlist(vc.traits.this.setting[which((grepl("Founder.pop", names(vc.traits.this.setting)))&
                                                                                 (grepl(paste("Rep.", this.rep,sep = ""), names(vc.traits.this.setting))))]),nrow = 4)
            this.genetic.value.var.covar <- matrix(unlist(vc.genetic.values.this.setting[which((grepl("Founder.pop", names(vc.traits.this.setting)))&
                                                                                                 (grepl(paste("Rep.", this.rep,sep = ""), names(vc.traits.this.setting))))]),nrow = 4)
            this.breeding.value.var.covar <- matrix(unlist(vc.breeding.values.this.setting[which((grepl("Founder.pop", names(vc.traits.this.setting)))&
                                                                                                   (grepl(paste("Rep.", this.rep,sep = ""), names(vc.traits.this.setting))))]),nrow = 4)
            
            trait.1.variance <- c(trait.1.variance, this.trait.var.covar[1,1])
            trait.2.variance <- c(trait.2.variance, this.trait.var.covar[2,2])
            trait.3.variance <- c(trait.3.variance, this.trait.var.covar[3,3])
            trait.4.variance <- c(trait.4.variance, this.trait.var.covar[4,4]) 
            genetic.value.1.variance <- c(genetic.value.1.variance, this.genetic.value.var.covar[1,1])
            genetic.value.2.variance <- c(genetic.value.2.variance, this.genetic.value.var.covar[2,2])
            genetic.value.3.variance <- c(genetic.value.3.variance, this.genetic.value.var.covar[3,3])
            genetic.value.4.variance <- c(genetic.value.4.variance, this.genetic.value.var.covar[4,4])
            breeding.value.1.variance <- c(breeding.value.1.variance, this.breeding.value.var.covar[1,1])
            breeding.value.2.variance <- c(breeding.value.2.variance, this.breeding.value.var.covar[2,2])
            breeding.value.3.variance <- c(breeding.value.3.variance, this.breeding.value.var.covar[3,3])
            breeding.value.4.variance <- c(breeding.value.4.variance, this.breeding.value.var.covar[4,4])        
            
            
            #For loop through the different types of selection
            
            
            for(this.selection.type in c("Direct","Disruptive","Stabilizing")){           
              for(this.selection.number in c("112", "224")){
                #vc.traits.this.subpop <- vc.traits.this.setting[which(grepl(this.selection.type, names(vc.traits.this.setting)))]
                vc.traits.this.subpop <- vc.traits.this.setting[which((grepl(this.selection.type, names(vc.traits.this.setting))) & 
                                                                        (grepl(this.selection.number, names(vc.traits.this.setting))) )]
                vc.genetic.values.this.subpop <- vc.genetic.values.this.setting[which((grepl(this.selection.type, names(vc.traits.this.setting))) & 
                                                                                        (grepl(this.selection.number, names(vc.traits.this.setting))) )]  
                vc.breeding.values.this.subpop <- vc.breeding.values.this.setting[which((grepl(this.selection.type, names(vc.traits.this.setting))) & 
                                                                                          (grepl(this.selection.number, names(vc.traits.this.setting))) )]  
                
                
                for(this.generation in 1:10){
                  
                  factor.A.vector <- c(factor.A.vector, factor.A)
                  factor.B.vector <- c(factor.B.vector, factor.B)
                  factor.C.vector <- c(factor.C.vector, factor.C)
                  factor.D.vector <- c(factor.D.vector, factor.D)
                  rep.vector <- c(rep.vector, this.rep)  
                  subpopulation.vector <- c(subpopulation.vector, this.selection.type)
                  selection.number.vector <- c(selection.number.vector, this.selection.number)
                  generation.vector <- c(generation.vector, this.generation)
                  
                  this.trait.var.covar <- matrix(unlist(vc.traits.this.subpop[which(grepl(paste("Gen.",this.generation, ".Factor", sep = ""), names(vc.traits.this.subpop)))]),nrow = 4)
                  this.genetic.value.var.covar <- matrix(unlist(vc.genetic.values.this.subpop[which(grepl(paste("Gen.",this.generation, ".Factor", sep = ""), names(vc.genetic.values.this.subpop)))]),nrow = 4)
                  this.breeding.value.var.covar <- matrix(unlist(vc.breeding.values.this.subpop[which(grepl(paste("Gen.",this.generation, ".Factor", sep = ""),names(vc.breeding.values.this.subpop)))]),nrow = 4)
                  
                  trait.1.variance <- c(trait.1.variance, this.trait.var.covar[1,1])
                  trait.2.variance <- c(trait.2.variance, this.trait.var.covar[2,2])
                  trait.3.variance <- c(trait.3.variance, this.trait.var.covar[3,3])
                  trait.4.variance <- c(trait.4.variance, this.trait.var.covar[4,4]) 
                  genetic.value.1.variance <- c(genetic.value.1.variance, this.genetic.value.var.covar[1,1])
                  genetic.value.2.variance <- c(genetic.value.2.variance, this.genetic.value.var.covar[2,2])
                  genetic.value.3.variance <- c(genetic.value.3.variance, this.genetic.value.var.covar[3,3])
                  genetic.value.4.variance <- c(genetic.value.4.variance, this.genetic.value.var.covar[4,4])
                  breeding.value.1.variance <- c(breeding.value.1.variance, this.breeding.value.var.covar[1,1])
                  breeding.value.2.variance <- c(breeding.value.2.variance, this.breeding.value.var.covar[2,2])
                  breeding.value.3.variance <- c(breeding.value.3.variance, this.breeding.value.var.covar[3,3])
                  breeding.value.4.variance <- c(breeding.value.4.variance, this.breeding.value.var.covar[4,4])
                  
                  
                }#End for loop through the generations:for(this.generation in 1:10)
              }#End for loop through the number of individuals to select
          }#End for loop through the different types of selection: for(this.selection type in c("Direct","Disruptive","Stabilizing"))
        }#End for(rep in 1:3)
      }#End for(el in 1:4)
    }#End for(k in 1:4)
  }#End for(j in 1:4)
}#End for(i in 1:3)

data.for.vc.scatterplot <- data.frame(trait.1.variance, trait.2.variance,
                                      trait.3.variance, trait.4.variance,
                                      genetic.value.1.variance, genetic.value.2.variance,
                                      genetic.value.3.variance, genetic.value.4.variance,
                                      breeding.value.1.variance, breeding.value.2.variance,
                                      breeding.value.3.variance, breeding.value.4.variance,
                                      factor.A.vector, factor.B.vector, factor.C.vector,
                                      factor.D.vector, rep.vector,subpopulation.vector, selection.number.vector, 
                                      generation.vector)



#For loop through the factors; index on "this.factor"
plot.trait.1 <- FALSE
plot.trait.2 <- TRUE
plot.trait.3 <- TRUE
plot.trait.4 <- TRUE
for(this.factor in which(grepl("factor",colnames(data.for.vc.scatterplot)))){
  #For loop through each level of the ith factor - there will be one pdf per factor (maybe create a subdirectory for this); index on "this.level"
    print(paste("this.factor = ",this.factor, sep = ""))  
  #For loop through the different levels of Factor i (columns)
     for(this.level in unique(data.for.vc.scatterplot[,this.factor])){
       data.for.vc.scatterplot.this.factor.level <- data.for.vc.scatterplot[which(data.for.vc.scatterplot[,this.factor] == this.level),]
       
       #################################Here is where the code for making the plot begins
 
       this.min.y.axis <- log(min(data.for.vc.scatterplot.this.factor.level[,which(grepl("variance",
                                                                                         colnames(data.for.vc.scatterplot.this.factor.level)))]))
       this.max.y.axis <- log(max(data.for.vc.scatterplot.this.factor.level[,which(grepl("variance",
                                                                                         colnames(data.for.vc.scatterplot.this.factor.level)))]))
      

       the.plot.trait.vc <- ggplot(data.for.vc.scatterplot.this.factor.level,  aes(x = generation.vector)) + ylim(this.min.y.axis,this.max.y.axis) + 
         xlab("Generation Number") + theme(axis.title.x = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20))+
         theme(axis.text = 8)

       the.plot.genetic.value.vc <- ggplot(data.for.vc.scatterplot.this.factor.level,  aes(x = generation.vector)) + ylim(this.min.y.axis,this.max.y.axis) + 
         xlab("Generation Number") + theme(axis.title.x = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20))+
         theme(axis.text = 8)

       the.plot.breeding.value.vc <- ggplot(data.for.vc.scatterplot.this.factor.level,  aes(x = generation.vector)) + ylim(this.min.y.axis,this.max.y.axis) + 
         xlab("Generation Number") + theme(axis.title.x = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20))+
         theme(axis.text = 8)
       
       line.types <- c("twodash", "solid", "longdash", "dotted", "dotdash", "dashed")
       count <- 1
       this.population <- NULL
       #For loop through the different kind of selection/subpopulation levels; index on this.selection.type 
       for(this.selection.type in unique(data.for.vc.scatterplot.this.factor.level$subpopulation.vector)[-1]){
         for(this.selection.number in unique(data.for.vc.scatterplot.this.factor.level$selection.number.vector)[-1]){
           this.line.type <- line.types[count]
           this.population <- c(this.population, paste(this.selection.type, ": ", this.selection.number, sep= ""))
           #Source in the code below that will help make these plots
           input.scatter.plot.data <- data.for.vc.scatterplot.this.factor.level[which((data.for.vc.scatterplot.this.factor.level$subpopulation.vector == "Founder")|
                                                                                        ((data.for.vc.scatterplot.this.factor.level$subpopulation.vector == this.selection.type)&
                                                                                           data.for.vc.scatterplot.this.factor.level$selection.number.vector == this.selection.number) ),] 
           
           #Add information for this subpopulation to the plot of variance components for traits
           source("/Users/alipka/Library/CloudStorage/Box-Box/Sabbatical_Roslin_Institute/R_workspace/omnigenic_sim/Functions_to_Make_Life_Easier/Make_VC_Scatter_Plots_Traits.R")
 
           #Add information for this subpopulation to the plot of variance components for the genetic values
           source("/Users/alipka/Library/CloudStorage/Box-Box/Sabbatical_Roslin_Institute/R_workspace/omnigenic_sim/Functions_to_Make_Life_Easier/Make_VC_Scatter_Plots_Genetic_Values.R")
           
           #Add information for this subpopulation to the plot of variance components for the breeding values
           source("/Users/alipka/Library/CloudStorage/Box-Box/Sabbatical_Roslin_Institute/R_workspace/omnigenic_sim/Functions_to_Make_Life_Easier/Make_VC_Scatter_Plots_Breeding_Values.R")
           
           count <- count+1
         }#End  for(this.selection.number in unique(data.for.vc.scatterplot.this.factor.level$selection.number.vector)[-1])
       }#End for(this.selection.type in unique(data.for.vc.scatterplot.this.factor.level$subpopulation.vector)[-1]) 
       

       pdf(paste("VC_Plots/VC.Factor",substr(colnames(data.for.vc.scatterplot)[this.factor],start = 7,stop = 9),"Level.",this.level,".plot.trait.VCs.pdf", sep = ""), width = 10)
       print(the.plot.trait.vc) 
       #This code is modified from here: https://stackoverflow.com/questions/57984046/how-to-create-a-stand-alone-legend-in-r
       plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
       legend("topleft", legend = this.population, lty = line.types)
       legend("topright", legend = c("Trait", "Core", "Peripheral", "CxP"),
              fill = c("black", "red", "blue", "purple"))
       dev.off()
       
       pdf(paste("VC_Plots/VC.Factor",substr(colnames(data.for.vc.scatterplot)[this.factor],start = 7,stop = 9),"Level.",this.level,".plot.genetic.value.VCs.pdf", sep = ""), width = 10)
       print(the.plot.genetic.value.vc)
       #This code is modified from here: https://stackoverflow.com/questions/57984046/how-to-create-a-stand-alone-legend-in-r
       plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
       legend("topleft", legend = this.population, lty = line.types)
       legend("topright", legend = c("Trait", "Core", "Peripheral", "CxP"),
              fill = c("black", "red", "blue", "purple"))
       dev.off()
       
       pdf(paste("VC_Plots/VC.Factor",substr(colnames(data.for.vc.scatterplot)[this.factor],start = 7,stop = 9),"Level.",this.level,".plot.breeding.value.VCs.pdf", sep = ""), width = 10)
       print(the.plot.breeding.value.vc)
       #This code is modified from here: https://stackoverflow.com/questions/57984046/how-to-create-a-stand-alone-legend-in-r
       plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
       legend("topleft", legend = this.population, lty = line.types)
       legend("topright", legend = c("Trait", "Core", "Peripheral", "CxP"),
              fill = c("black", "red", "blue", "purple"))
       dev.off()
    }#End for(this.level in unique(data.for.vc.scatterplot[,this.factor]))for loop through each level of the ith factor
}#End  for(this.factor in grepl("factor",colnames(data.for.vc.scatterplot))) for loop through the factors



#######Create a very similar set of plots; this time with only one combination of all four factors
plot.trait.1 <- FALSE
plot.trait.2 <- TRUE
plot.trait.3 <- TRUE
plot.trait.4 <- TRUE

this.level.of.factor.A <- 1
this.level.of.factor.B <- 0.05
this.level.of.factor.C <- 0.05
this.level.of.factor.D <- 0.05

for(this.level.of.factor.A in unique(data.for.vc.scatterplot[,13])){
  for(this.level.of.factor.B in unique(data.for.vc.scatterplot[,14])){
    for(this.level.of.factor.C in unique(data.for.vc.scatterplot[,15])){
      for(this.level.of.factor.D in unique(data.for.vc.scatterplot[,16])){
  
        data.for.vc.scatterplot.this.factor.level <- data.for.vc.scatterplot[which((data.for.vc.scatterplot[,13] == this.level.of.factor.A)&
                                                                                   (data.for.vc.scatterplot[,14] == this.level.of.factor.B)&
                                                                                   (data.for.vc.scatterplot[,15] == this.level.of.factor.C)&
                                                                                   (data.for.vc.scatterplot[,16] == this.level.of.factor.D)),]
        
        #################################Here is where the code for making the plot begins
        
        this.min.y.axis <- log(min(data.for.vc.scatterplot.this.factor.level[,which(grepl("variance",
                                                                                          colnames(data.for.vc.scatterplot.this.factor.level)))]))
        this.max.y.axis <- log(max(data.for.vc.scatterplot.this.factor.level[,which(grepl("variance",
                                                                                          colnames(data.for.vc.scatterplot.this.factor.level)))]))
        
        
        the.plot.trait.vc <- ggplot(data.for.vc.scatterplot.this.factor.level,  aes(x = generation.vector)) + ylim(this.min.y.axis,this.max.y.axis) + 
          xlab("Generation Number") + theme(axis.title.x = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20))+
          theme(axis.text = element_text(size = 15))
        
        the.plot.genetic.value.vc <- ggplot(data.for.vc.scatterplot.this.factor.level,  aes(x = generation.vector)) + ylim(this.min.y.axis,this.max.y.axis) + 
          xlab("Generation Number") + theme(axis.title.x = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20))+
          theme(axis.text = element_text(size = 15))
        
        the.plot.breeding.value.vc <- ggplot(data.for.vc.scatterplot.this.factor.level,  aes(x = generation.vector)) + ylim(this.min.y.axis,this.max.y.axis) + 
          xlab("Generation Number") + theme(axis.title.x = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20))+
          theme(axis.text = element_text(size = 15))
        
        line.types <- c("twodash", "solid", "longdash", "dotted", "dotdash", "dashed")
        count <- 1
        this.population <- NULL
        #For loop through the different kind of selection/subpopulation levels; index on this.selection.type 
        for(this.selection.type in unique(data.for.vc.scatterplot.this.factor.level$subpopulation.vector)[-1]){
          for(this.selection.number in unique(data.for.vc.scatterplot.this.factor.level$selection.number.vector)[-1]){
            this.line.type <- line.types[count]
            this.population <- c(this.population, paste(this.selection.type, ": ", this.selection.number, sep= ""))
            #Source in the code below that will help make these plots
            input.scatter.plot.data <- data.for.vc.scatterplot.this.factor.level[which((data.for.vc.scatterplot.this.factor.level$subpopulation.vector == "Founder")|
                                                                                         ((data.for.vc.scatterplot.this.factor.level$subpopulation.vector == this.selection.type)&
                                                                                            data.for.vc.scatterplot.this.factor.level$selection.number.vector == this.selection.number) ),] 
            
            #Add information for this subpopulation to the plot of variance components for traits
            source("/Users/alipka/Library/CloudStorage/Box-Box/Sabbatical_Roslin_Institute/R_workspace/omnigenic_sim/Functions_to_Make_Life_Easier/Make_VC_Scatter_Plots_Traits.R")
            
            #Add information for this subpopulation to the plot of variance components for the genetic values
            source("/Users/alipka/Library/CloudStorage/Box-Box/Sabbatical_Roslin_Institute/R_workspace/omnigenic_sim/Functions_to_Make_Life_Easier/Make_VC_Scatter_Plots_Genetic_Values.R")
            
            #Add information for this subpopulation to the plot of variance components for the breeding values
            source("/Users/alipka/Library/CloudStorage/Box-Box/Sabbatical_Roslin_Institute/R_workspace/omnigenic_sim/Functions_to_Make_Life_Easier/Make_VC_Scatter_Plots_Breeding_Values.R")
            
            count <- count+1
          }#End  for(this.selection.number in unique(data.for.vc.scatterplot.this.factor.level$selection.number.vector)[-1])
        }#End for(this.selection.type in unique(data.for.vc.scatterplot.this.factor.level$subpopulation.vector)[-1]) 
        
        
        pdf(paste("VC_Plots/VC_Plots_by_Trt_Level/VC.", paste("Factor.A.", this.level.of.factor.A,
                                                              ".Factor.B.", this.level.of.factor.B,
                                                              "Factor.C.", this.level.of.factor.C,
                                                              "Factor.D.", this.level.of.factor.D,
                                                              sep = ""), ".plot.trait.VCs.pdf", sep = ""), width = 10)
        print(the.plot.trait.vc) 
        #This code is modified from here: https://stackoverflow.com/questions/57984046/how-to-create-a-stand-alone-legend-in-r
        plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
        legend("topleft", legend = this.population, lty = line.types)
        legend("topright", legend = c("Trait", "Core", "Peripheral", "CxP"),
               fill = c("black", "red", "blue", "purple"))
        dev.off()
        
        pdf(paste("VC_Plots/VC_Plots_by_Trt_Level/VC.", paste("Factor.A.", this.level.of.factor.A,
                                                              ".Factor.B.", this.level.of.factor.B,
                                                              "Factor.C.", this.level.of.factor.C,
                                                              "Factor.D.", this.level.of.factor.D,
                                                              sep = ""), ".plot.genetic.value.VCs.pdf", sep = ""), width = 10)
        print(the.plot.genetic.value.vc)
        #This code is modified from here: https://stackoverflow.com/questions/57984046/how-to-create-a-stand-alone-legend-in-r
        plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
        legend("topleft", legend = this.population, lty = line.types)
        legend("topright", legend = c("Trait", "Core", "Peripheral", "CxP"),
               fill = c("black", "red", "blue", "purple"))
        dev.off()
        
        pdf(paste("VC_Plots/VC_Plots_by_Trt_Level/VC.", paste("Factor.A.", this.level.of.factor.A,
                                                             ".Factor.B.", this.level.of.factor.B,
                                                             "Factor.C.", this.level.of.factor.C,
                                                             "Factor.D.", this.level.of.factor.D,
                                                             sep = ""),".plot.breeding.value.VCs.pdf", sep = ""), width = 10)
        print(the.plot.breeding.value.vc)
        #This code is modified from here: https://stackoverflow.com/questions/57984046/how-to-create-a-stand-alone-legend-in-r
        plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
        legend("topleft", legend = this.population, lty = line.types)
        legend("topright", legend = c("Trait", "Core", "Peripheral", "CxP"),
               fill = c("black", "red", "blue", "purple"))
        dev.off()
        
      }#end for(this.level.of.factor.D in unique(data.for.vc.scatterplot[,16]))
    } #end unique(data.for.vc.scatterplot[,15])
  } # end for(this.level.of.factor.B in unique(data.for.vc.scatterplot[,14]))
}#end for(this.level.of.factor.A in unique(data.for.vc.scatterplot[,13]))




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


##############Same analysis, using SNPs instead of QTLs
these.spearman.rank.correlation.between.GWAS.peripheral.SNPs <- readRDS("/Users/alipka/Library/CloudStorage/Box-Box/IR-281/Master-Result-Lists-2024-06-13/master.these.spearman.rank.correlation.between.GWAS.peripheral.SNPs.RDS")

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
          results.this.setting <- these.spearman.rank.correlation.between.GWAS.peripheral.SNPs[which(names(these.spearman.rank.correlation.between.GWAS.peripheral.QTNs)
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
pdf("GWAS.Results.Peripheral.SNPs.written.while.writing.plot.code.pdf", width = 8)
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
#####Using QTLs
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
pdf("Box.plot.of.Prediciton.Accuracies.using.QTNs.20240710.pdf", width = 50)
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
    source("/Users/alipka/Library/CloudStorage/Box-Box/Sabbatical_Roslin_Institute/R_workspace/Sabbatical_Project/Functions_to_Make_Life_Easier/Make_GS_Box_Plots.R")
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
    source("/Users/alipka/Library/CloudStorage/Box-Box/Sabbatical_Roslin_Institute/R_workspace/Sabbatical_Project/Functions_to_Make_Life_Easier/Make_GS_Box_Plots.R")
    fa.row.2 <- fa
    fb.row.2 <- fb
    fc.row.2 <- fc
    fd.row.2 <- fd
    
    #########Make the third row of the plot
    #Source in the code below that will help make these plots
    input.box.plot.data <- data.for.boxplot.this.validation.pop.melted.both.pop.train
    y.axis.label <- paste("Pop ",these.two.training.sets[1], "+", these.two.training.sets[2], " TS", sep = "") 
    
    #Source in some gplot code that will make the plots, and save them as objects
    source("/Users/alipka/Library/CloudStorage/Box-Box/Sabbatical_Roslin_Institute/R_workspace/Sabbatical_Project/Functions_to_Make_Life_Easier/Make_GS_Box_Plots.R")
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


#####Using SNPs
# maybe consider putting this into a function
these.prediction.accuracies.SNPs <- readRDS("/Users/alipka/Library/CloudStorage/Box-Box/IR-281/Master-Result-Lists-2024-06-13/master.these.prediction.accuracies.SNPs.RDS")

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
          results.this.setting <- these.prediction.accuracies.SNPs[which(names(these.prediction.accuracies.SNPs)
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
pdf("Box.plot.of.Prediciton.Accuracies.using.SNPs.20240710.pdf", width = 50)
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
  source("/Users/alipka/Library/CloudStorage/Box-Box/Sabbatical_Roslin_Institute/R_workspace/Sabbatical_Project/Functions_to_Make_Life_Easier/Make_GS_Box_Plots.R")
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
  source("/Users/alipka/Library/CloudStorage/Box-Box/Sabbatical_Roslin_Institute/R_workspace/Sabbatical_Project/Functions_to_Make_Life_Easier/Make_GS_Box_Plots.R")
  fa.row.2 <- fa
  fb.row.2 <- fb
  fc.row.2 <- fc
  fd.row.2 <- fd
  
  #########Make the third row of the plot
  #Source in the code below that will help make these plots
  input.box.plot.data <- data.for.boxplot.this.validation.pop.melted.both.pop.train
  y.axis.label <- paste("Pop ",these.two.training.sets[1], "+", these.two.training.sets[2], " TS", sep = "") 
  
  #Source in some gplot code that will make the plots, and save them as objects
  source("/Users/alipka/Library/CloudStorage/Box-Box/Sabbatical_Roslin_Institute/R_workspace/Sabbatical_Project/Functions_to_Make_Life_Easier/Make_GS_Box_Plots.R")
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

