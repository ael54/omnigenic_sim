


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
#Read in the Spearman correlation coefficients between GWAS results
#correlation.core.QTNs <- readRDS("/Users/alipka/Library/CloudStorage/Box-Box/IR-281/results-analysis-20250416/master.these.spearman.rank.correlation.between.GWAS.core.QTNs.RDS") 
#correlation.core.GWAS <- readRDS("/Users/alipka/Library/CloudStorage/Box-Box/IR-281/results-analysis-20250416/master.these.spearman.rank.correlation.between.GWAS.core.SNPs.RDS")
#correlation.peripheral.QTNs <- readRDS("/Users/alipka/Library/CloudStorage/Box-Box/IR-281/results-analysis-20250416/master.these.spearman.rank.correlation.between.GWAS.peripheral.QTNs.RDS") 
#correlation.peripheral.GWAS <- readRDS("/Users/alipka/Library/CloudStorage/Box-Box/IR-281/results-analysis-20250416/master.these.spearman.rank.correlation.between.GWAS.peripheral.SNPs.RDS")

correlation.core.QTNs <- readRDS("/Users/alipka/Library/CloudStorage/Box-Box/IR-281/results-analysis-20250416/master.these.spearman.rank.correlation.between.GWAS.core.QTNs.RDS") 
correlation.core.GWAS <- readRDS("/Users/alipka/Library/CloudStorage/Box-Box/IR-281/results-analysis-20250416/master.these.spearman.rank.correlation.between.GWAS.core.SNPs.RDS")
correlation.peripheral.QTNs <- readRDS("/Users/alipka/Library/CloudStorage/Box-Box/IR-281/results-analysis-20250416/master.these.spearman.rank.correlation.between.GWAS.peripheral.QTNs.RDS") 
correlation.peripheral.GWAS <- readRDS("/Users/alipka/Library/CloudStorage/Box-Box/IR-281/results-analysis-20250416/master.these.spearman.rank.correlation.between.GWAS.peripheral.SNPs.RDS")


#Format the results for trait values, genetic values, and breeding values so that they can be
# used to make graphs



#For each simulation setting: extract the Spearman rank correlation coefficients
# for core QTNs, core SNPs, peripheral QTNs, and peripheral SNPs
QTN.or.SNP <- NULL
core.or.peripheral <- NULL 
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
            #This code will work once all of the simulation settings are populated
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
            #END: This code will work once all of the simulation settings are populated
            
            for(sim.rep.index in 1:3){
              correlation.core.QTNs.this.setting <- correlation.core.QTNs[[sim.rep.index]][which(names(correlation.core.QTNs[[sim.rep.index]])!= "NA")]
              correlation.peripheral.QTNs.this.setting <- correlation.peripheral.QTNs[[sim.rep.index]][which(names(correlation.peripheral.QTNs[[sim.rep.index]])!= "NA")]
              if(names(correlation.core.QTNs.this.setting)!=names(correlation.peripheral.QTNs.this.setting)) next
            
              #Keep track of which setting we are looking at
              this.actual.setting <- names(correlation.core.QTNs.this.setting)
              this.setting <- this.actual.setting
              
              #Convert the core and peripheral results to a data frame so we can make graphs
              correlation.core.QTNs.this.setting <- as.data.frame(correlation.core.QTNs.this.setting)
              colnames(correlation.core.QTNs.this.setting) <- c("Pop1", "Pop2", "Spear.Cor.Core")
              
              correlation.peripheral.QTNs.this.setting <- as.data.frame(correlation.peripheral.QTNs.this.setting)
              colnames(correlation.peripheral.QTNs.this.setting) <- c("Pop1", "Pop2", "Spear.Cor.Peri")
              
              #Loop through each of the six target subpopulations, and make a single graph for each of them
              pdf(paste("GWAS_Plots/Spearman.Cor.",this.setting,".20250422.pdf",sep = ""), width = 12)
              for(this.plot in unique(correlation.peripheral.QTNs.this.setting$Pop1)){
                if(grepl(paste(".prev.gen", sep = ""), this.plot)) next 
                
                core.correlations.this.subpop <- correlation.core.QTNs.this.setting[which((correlation.core.QTNs.this.setting$Pop1 == this.plot)|
                                                                                          (correlation.core.QTNs.this.setting$Pop2 == this.plot)),]
                peripheral.correlations.this.subpop <- correlation.peripheral.QTNs.this.setting[which((correlation.peripheral.QTNs.this.setting$Pop1 == this.plot)|
                                                                                                  (correlation.core.QTNs.this.setting$Pop2 == this.plot)),]
                
                #Relabel the previous generation to the name of the subpopulation under study
                core.correlations.this.subpop$Pop2[which(grepl(this.plot, core.correlations.this.subpop$Pop2))] = this.plot
                peripheral.correlations.this.subpop$Pop2[which(grepl(this.plot, peripheral.correlations.this.subpop$Pop2))] = this.plot
                
                #For all but the first setting, There will be some targeted populations that are labeled in Pop 1 due to the naming scheme
                # Move labels over to Pop2
                if(length(core.correlations.this.subpop$Pop1 != this.plot) >0 ){
                  core.correlations.this.subpop$Pop2[which(core.correlations.this.subpop$Pop1 != this.plot)] =
                    core.correlations.this.subpop$Pop1[which(core.correlations.this.subpop$Pop1 != this.plot)]
                  peripheral.correlations.this.subpop$Pop2[which(peripheral.correlations.this.subpop$Pop1 != this.plot)] =
                    peripheral.correlations.this.subpop$Pop1[which(peripheral.correlations.this.subpop$Pop1 != this.plot)]
                }# end if(length(core.correlations.this.subpop$Pop1 != this plot) >0 )
                
                #Get rid of all "prev.gen" levels in Pop.2
                core.correlations.this.subpop <- core.correlations.this.subpop[-which(grepl(paste(".prev.gen"), core.correlations.this.subpop$Pop2)),]
                peripheral.correlations.this.subpop <- peripheral.correlations.this.subpop[-which(grepl(paste(".prev.gen"), peripheral.correlations.this.subpop$Pop2)),]
                
                #Set Pop2 to a factor
                core.correlations.this.subpop$Pop2 <- factor(core.correlations.this.subpop$Pop2, levels = c("Directional.selection.10.pct", 
                                                                                                               "Directional.selection.20.pct", 
                                                                                                               "Disruptive.selection.10.pct", 
                                                                                                               "Disruptive.selection.20.pct", 
                                                                                                               "Stabilizing.selection.10.pct", 
                                                                                                               "Stabilizing.selection.20.pct"))
                peripheral.correlations.this.subpop$Pop2 <- factor(peripheral.correlations.this.subpop$Pop2, levels = c("Directional.selection.10.pct", 
                                                                                                            "Directional.selection.20.pct", 
                                                                                                            "Disruptive.selection.10.pct", 
                                                                                                            "Disruptive.selection.20.pct", 
                                                                                                            "Stabilizing.selection.10.pct", 
                                                                                                            "Stabilizing.selection.20.pct"))
                
                #Make the plots. Turn the code below into a separate R function

                
                the.merged.data.for.plot <- merge(core.correlations.this.subpop, peripheral.correlations.this.subpop, by.x = "Pop2", by.y = "Pop2")
                
                this.min.y.axis <- -1
                this.max.y.axis <- 1
                #Initiate the plot
                the.plot.QTN.correlations <-   ggplot(the.merged.data.for.plot,  aes(x = Pop2)) + ylim(this.min.y.axis,this.max.y.axis) +
                                      xlab("Subpopulation") +theme(axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5))+
                                      theme(theme(axis.title.x = element_text(size = 20)))+                  
                                      ylab("Spearman Correlation") + theme(axis.title.y = element_text(size = 20))+
                                      theme(axis.title.x = element_text(size = 20)) + ggtitle(paste(this.setting, ": ", this.plot, sep = ""))
                                                        
                
                #Add information on the core genes
                the.plot.QTN.correlations <- the.plot.QTN.correlations + geom_point(data= the.merged.data.for.plot, aes(x = Pop2,y=  Spear.Cor.Core), size = 3,  colour = "#CC0000" )
                
                #Add information on the peripheral genes
                the.plot.QTN.correlations <- the.plot.QTN.correlations + geom_point(data= the.merged.data.for.plot, aes(x = Pop2, y=  Spear.Cor.Peri), size = 3,  colour = "#0000FF" )
                
                #geom_point(col = "blue") + geom_point(data= the.merged.data.for.plot, aes(x = Pop2, y= Spear.Cor.Peri, col ="#6600FF"))+
                #   + +
                  
                                 #theme(axis.title.x = element_text(size = 20))
                  #xlab("Subpopulation") + theme(axis.title.x = element_text(size = 10)) + theme(axis.title.y = element_text(size = 20),)+
                  #theme(axis.text = 8)
                
                plot(the.plot.QTN.correlations)
                
        
              
                
              }#end for(this.plot in unique(correlation.peripheral.QTNs.this.setting$Pop1))
              dev.off()
              
                
             
              
            }# end for(sim.rep.index ...)
            #old code for emulating
            vc.traits.this.setting <- vc.traits[which(grepl(this.setting, names(vc.traits)))]
            vc.genetic.values.this.setting <- vc.genetic.values[which(grepl(this.setting, names(vc.genetic.values)))]  
            vc.breeding.values.this.setting <- vc.breeding.values[which(grepl(this.setting, names(vc.breeding.values)))]  
            #end old code for emulating
            
            
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



