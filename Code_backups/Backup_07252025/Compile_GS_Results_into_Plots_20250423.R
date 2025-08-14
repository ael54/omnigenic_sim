


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

prediction.accuracies.QTNs <- readRDS("/Users/alipka/Library/CloudStorage/Box-Box/IR-281/results-analysis-20250421/master.these.prediction.accuracies.QTNs.RDS") 
prediction.accuracies.SNPs <- readRDS("/Users/alipka/Library/CloudStorage/Box-Box/IR-281/results-analysis-20250421/master.these.prediction.accuracies.SNPs.RDS")


#Format the results for trait values, genetic values, and breeding values so that they can be
# used to make graphs



#For each simulation setting: extract the Spearman rank correlation coefficients
# for core QTNs, core SNPs, peripheral QTNs, and peripheral SNPs
     
  for(sim.rep.index in 1:length(correlation.core.QTNs)){
    
    this.setting <- names(correlation.core.QTNs[sim.rep.index][1][[1]])
    correlation.core.QTNs.this.setting <- correlation.core.QTNs[sim.rep.index][1][[1]]
    names(correlation.core.QTNs.this.setting) <- paste("Core.Results")
    this.core.pop1 <- correlation.core.QTNs.this.setting$Core.Results$the.pop.1
    this.core.pop2 <- correlation.core.QTNs.this.setting$Core.Results$the.pop.2
    this.core.spearman <- correlation.core.QTNs.this.setting$Core.Results$the.spearman.rank
    correlation.core.QTNs.this.setting <- data.frame(this.core.pop1,this.core.pop2,
                                                     this.core.spearman)
    
    colnames(correlation.core.QTNs.this.setting) <- c("Pop1", "Pop2", "Spear.Cor.Core")
    
    
    this.setting.peripheral <- names(correlation.peripheral.QTNs[sim.rep.index][1][[1]])
    if(this.setting != this.setting.peripheral) next
    
    correlation.peripheral.QTNs.this.setting <- correlation.peripheral.QTNs[sim.rep.index][1][[1]]
    names(correlation.peripheral.QTNs.this.setting) <- paste("peripheral.Results")
    this.peripheral.pop1 <- correlation.peripheral.QTNs.this.setting$peripheral.Results$the.pop.1
    this.peripheral.pop2 <- correlation.peripheral.QTNs.this.setting$peripheral.Results$the.pop.2
    this.peripheral.spearman <- correlation.peripheral.QTNs.this.setting$peripheral.Results$the.spearman.rank
    correlation.peripheral.QTNs.this.setting <- data.frame(this.peripheral.pop1,this.peripheral.pop2,
                                                           this.peripheral.spearman)
    
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




