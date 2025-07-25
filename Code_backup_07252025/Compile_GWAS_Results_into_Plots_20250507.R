


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
library(tidyverse)



#########

correlation.core.QTNs <- readRDS("/Users/alipka/Library/CloudStorage/Box-Box/IR-281/results-analysis-20250421/master.these.spearman.rank.correlation.between.GWAS.core.QTNs.RDS") 
correlation.core.SNPs <- readRDS("/Users/alipka/Library/CloudStorage/Box-Box/IR-281/results-analysis-20250421/master.these.spearman.rank.correlation.between.GWAS.core.SNPs.RDS")
correlation.peripheral.QTNs <- readRDS("/Users/alipka/Library/CloudStorage/Box-Box/IR-281/results-analysis-20250421/master.these.spearman.rank.correlation.between.GWAS.peripheral.QTNs.RDS") 
correlation.peripheral.SNPs <- readRDS("/Users/alipka/Library/CloudStorage/Box-Box/IR-281/results-analysis-20250421/master.these.spearman.rank.correlation.between.GWAS.peripheral.SNPs.RDS")


#Format the results for trait values, genetic values, and breeding values so that they can be
# used to make graphs



#For each simulation setting: extract the Spearman rank correlation coefficients
# for core QTNs, core SNPs, peripheral QTNs, and peripheral SNPs
  sim.rep.index <- 1 
  is.rep.1 <-  TRUE
  #in 1:length(correlation.core.QTNs)
  while(is.rep.1 ==TRUE){
    
    
    this.setting.rep.1 <- names(correlation.core.QTNs[sim.rep.index][1][[1]])
    #If we are beyond rep 1, then stop this while loop
    this.rep.number <- as.numeric(substr(this.setting.rep.1, start = nchar(this.setting.rep.1) - 4, 
                              stop = nchar(this.setting.rep.1) - 4))
    if(this.rep.number != 1){ 
      print(paste("Gee whiz! This while loop stopped on ", sim.rep.index, sep = ""))
      is.rep.1 == FALSE
      stop
    }#End if(this.rep.number != 1)
    
    correlation.core.QTNs.this.setting.rep.1 <- correlation.core.QTNs[sim.rep.index][1][[1]]
    names(correlation.core.QTNs.this.setting.rep.1) <- paste("Core.Results")
    this.core.pop1 <- correlation.core.QTNs.this.setting.rep.1$Core.Results$the.pop.1
    this.core.pop2 <- correlation.core.QTNs.this.setting.rep.1$Core.Results$the.pop.2
    this.core.spearman <- correlation.core.QTNs.this.setting.rep.1$Core.Results$the.spearman.rank
    correlation.core.QTNs.this.setting.rep.1 <- data.frame(this.core.pop1,this.core.pop2,
                                                     this.core.spearman)
    
    colnames(correlation.core.QTNs.this.setting.rep.1) <- c("Pop1", "Pop2", "Spear.Cor.Core")
    
    
    this.setting.peripheral.rep.1 <- names(correlation.peripheral.QTNs[sim.rep.index][1][[1]])
    if(this.setting.rep.1 != this.setting.peripheral.rep.1) next
    
    correlation.peripheral.QTNs.this.setting.rep.1 <- correlation.peripheral.QTNs[sim.rep.index][1][[1]]
    names(correlation.peripheral.QTNs.this.setting.rep.1) <- paste("peripheral.Results")
    this.peripheral.pop1 <- correlation.peripheral.QTNs.this.setting.rep.1$peripheral.Results$the.pop.1
    this.peripheral.pop2 <- correlation.peripheral.QTNs.this.setting.rep.1$peripheral.Results$the.pop.2
    this.peripheral.spearman <- correlation.peripheral.QTNs.this.setting.rep.1$peripheral.Results$the.spearman.rank
    correlation.peripheral.QTNs.this.setting.rep.1 <- data.frame(this.peripheral.pop1,this.peripheral.pop2,
                                                           this.peripheral.spearman)
    
    colnames(correlation.peripheral.QTNs.this.setting.rep.1) <- c("Pop1", "Pop2", "Spear.Cor.Peri")
    
    this.setting.factor.levels <- substr(this.setting.rep.1, start = 1 , stop = nchar(this.setting.rep.1) - 5)
    
    rep.index <- NULL
    for(this.iteration in (sim.rep.index+1):length(correlation.core.QTNs)){
      this.test.setting <- names(correlation.peripheral.QTNs[this.iteration][1][[1]])
      this.test.setting.factor.levels <- substr(this.test.setting, start = 1 , stop = nchar(this.test.setting) - 5)
      if(this.test.setting.factor.levels == this.setting.factor.levels){
        print(paste("The next rep is on iteration, ", this.iteration, sep = ""))
        rep.index <- c(rep.index, this.iteration) 
      } # End if(this.test.setting.factor.levels == this.setting.factor.levels)
    }#End for(this.iteration in 1:length(correlation.core.QTNs))
    
    correlation.core.QTNs.this.setting <- correlation.core.QTNs.this.setting.rep.1
    correlation.peripheral.QTNs.this.setting <- correlation.peripheral.QTNs.this.setting.rep.1
    for(other.reps in rep.index){
      #Create all of the information needed for core QTNs for the next reps
      this.setting.next.rep <- names(correlation.core.QTNs[other.reps][1][[1]])
      correlation.core.QTNs.this.setting.next.rep <- correlation.core.QTNs[other.reps][1][[1]]
      names(correlation.core.QTNs.this.setting.next.rep) <- paste("Core.Results")
      this.core.pop1 <- correlation.core.QTNs.this.setting.next.rep$Core.Results$the.pop.1
      this.core.pop2 <- correlation.core.QTNs.this.setting.next.rep$Core.Results$the.pop.2
      this.core.spearman <- correlation.core.QTNs.this.setting.next.rep$Core.Results$the.spearman.rank
      correlation.core.QTNs.this.setting.next.rep <- data.frame(this.core.pop1,this.core.pop2,
                                                             this.core.spearman)
      colnames(correlation.core.QTNs.this.setting.next.rep) <- c("Pop1", "Pop2", "Spear.Cor.Core")
      correlation.core.QTNs.this.setting <- rbind(correlation.core.QTNs.this.setting, correlation.core.QTNs.this.setting.next.rep)
      
      #Create all of the information needed for peripheral QTNs for the next reps
      this.setting.peripheral.next.rep <- names(correlation.peripheral.QTNs[other.reps][1][[1]])
      if(this.setting.next.rep != this.setting.peripheral.next.rep){
        next
        print(paste("Oh no!!!!!! The setting names did not match for simulation setting number ", other.reps, sep = ""))
      } 
      
      correlation.peripheral.QTNs.this.setting.next.rep <- correlation.peripheral.QTNs[other.reps][1][[1]]
      names(correlation.peripheral.QTNs.this.setting.next.rep) <- paste("peripheral.Results")
      this.peripheral.pop1 <- correlation.peripheral.QTNs.this.setting.next.rep$peripheral.Results$the.pop.1
      this.peripheral.pop2 <- correlation.peripheral.QTNs.this.setting.next.rep$peripheral.Results$the.pop.2
      this.peripheral.spearman <- correlation.peripheral.QTNs.this.setting.next.rep$peripheral.Results$the.spearman.rank
      correlation.peripheral.QTNs.this.setting.next.rep<- data.frame(this.peripheral.pop1,this.peripheral.pop2,
                                                                   this.peripheral.spearman)
      
      colnames(correlation.peripheral.QTNs.this.setting.next.rep) <- c("Pop1", "Pop2", "Spear.Cor.Peri")
      
      correlation.peripheral.QTNs.this.setting <- rbind(correlation.peripheral.QTNs.this.setting, correlation.peripheral.QTNs.this.setting.next.rep)
      
    }#End for(other.reps in rep.index)
    
    #Make the column names of correlations match up between core and peripheral genes
    colnames(correlation.core.QTNs.this.setting)[3] <- "Correlation"
    colnames(correlation.peripheral.QTNs.this.setting)[3] <- "Correlation"
    
    #Add a column in each of the core and peripheral genes files that indicate whether the 
    # QTNs are core or peripheral
    
    correlation.core.QTNs.this.setting <- cbind(correlation.core.QTNs.this.setting, 
                                                rep("Core", nrow(correlation.core.QTNs.this.setting)))
    colnames(correlation.core.QTNs.this.setting)[4] <- "Gene"

    correlation.peripheral.QTNs.this.setting <- cbind(correlation.peripheral.QTNs.this.setting, 
                                                rep("Peripheral", nrow(correlation.peripheral.QTNs.this.setting)))
    colnames(correlation.peripheral.QTNs.this.setting)[4] <- "Gene"
    
    core.and.peri.data <- rbind(correlation.core.QTNs.this.setting, correlation.peripheral.QTNs.this.setting)
    core.and.peri.data$Gene <- as.factor(core.and.peri.data$Gene)
    
    count <- 0
    for(this.comparison.population in unique(core.and.peri.data$Pop1)){
      if(grepl(paste(".prev.gen", sep = ""), this.comparison.population)) next 
      
      correlations.for.this.subpop <- core.and.peri.data[which((core.and.peri.data$Pop1 == this.comparison.population)|
                                                                                  (core.and.peri.data$Pop2 == this.comparison.population)),]
        
      #Relabel the previous generation to the name of the subpopulation under study
      correlations.for.this.subpop$Pop2[which(grepl(this.comparison.population, correlations.for.this.subpop$Pop2))] = this.comparison.population
       
      
      #For all but the first setting, There will be some targeted populations that are labeled in Pop 1 due to the naming scheme
      # Move labels over to Pop2
      if(length(which(correlations.for.this.subpop$Pop1 != this.comparison.population)) >0 ){
        #Swap the labels of pop 2 and pop 1
        these.indices <- which(correlations.for.this.subpop$Pop1 != this.comparison.population)
        move.these.Pop2.labels.to.Pop1 <- correlations.for.this.subpop$Pop2[these.indices]
        move.these.Pop1.labels.to.Pop2 <- correlations.for.this.subpop$Pop1[these.indices]
        correlations.for.this.subpop$Pop1[these.indices] <- move.these.Pop2.labels.to.Pop1
        correlations.for.this.subpop$Pop2[these.indices] <- move.these.Pop1.labels.to.Pop2
      }# end if(length(core.correlations.this.subpop$Pop1 != this plot) >0 )
      
      #Get rid of all "prev.gen" levels in Pop.2
      correlations.for.this.subpop <- correlations.for.this.subpop[-which(grepl(paste(".prev.gen"), correlations.for.this.subpop$Pop2)),]
       
      if(count == 0){
        input.data.for.boris.plot <- correlations.for.this.subpop
      }else{
        input.data.for.boris.plot <- rbind(input.data.for.boris.plot, correlations.for.this.subpop)
      }#end if(count == 0)
      
      count <- count + 1
    
    } #End for(this.comparison.population in unique(cor.and.peri.data$Pop1))
    
    #Make population and comparison population are factors
    colnames(input.data.for.boris.plot)[c(1,2)] <- c("Comp_Population", "Target_Population")
    input.data.for.boris.plot$Comp_Population <- factor(input.data.for.boris.plot$Comp_Population, levels = c("Directional.selection.10.pct", 
                                                          "Directional.selection.20.pct", 
                                                          "Disruptive.selection.10.pct", 
                                                          "Disruptive.selection.20.pct", 
                                                          "Stabilizing.selection.10.pct", 
                                                          "Stabilizing.selection.20.pct"))
    
    input.data.for.boris.plot$Target_Population <- factor(input.data.for.boris.plot$Target_Population, levels = c("Directional.selection.10.pct", 
                                                                                                              "Directional.selection.20.pct", 
                                                                                                              "Disruptive.selection.10.pct", 
                                                                                                              "Disruptive.selection.20.pct", 
                                                                                                              "Stabilizing.selection.10.pct", 
                                                                                                              "Stabilizing.selection.20.pct"))
    
    
    squared.correlation <- (input.data.for.boris.plot$Correlation)^2  
    input.data.for.boris.plot <- cbind(input.data.for.boris.plot, squared.correlation)
    #Make Boris' cool plot
    
    p1 <- ggplot(input.data.for.boris.plot, aes(x = Comp_Population, y = squared.correlation)) +
      geom_boxplot(aes(fill = Target_Population, linetype = Gene)) +
      labs(x = "Comparison population", y = "Spearman correlation squared") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    pdf(paste("GWAS_Plots/Box.plot.",this.setting.rep.1,".20250509.pdf",sep = ""), width = 12)
      print(p1)
    dev.off()
    
    if(sim.rep.index ==  1){
      data.for.overall.boris.plot <- input.data.for.boris.plot 
    } else{
      data.for.overall.boris.plot <- rbind(data.for.overall.boris.plot, input.data.for.boris.plot)
    }
    
    sim.rep.index <- sim.rep.index + 1
    rep.index <- NULL
  }# end for(sim.rep.index ...)  
    ####################################################Old code
 
  #Make an overall box plot
  p2 <- ggplot(data.for.overall.boris.plot, aes(x = Comp_Population, y = squared.correlation)) +
    geom_boxplot(aes(fill = Target_Population, linetype = Gene)) +
    labs(x = "Comparison population", y = "Spearman correlation squared") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  pdf(paste("GWAS_Plots/An.Overall.Box.plot.GWAS.Correlations.20250509.pdf",sep = ""), width = 12)
    print(p2)
  dev.off()
  
  
  #Make a plot like the last slide that Boris Mad
  p3 <- ggplot(data.for.overall.boris.plot, aes(x = Comp_Population, y = squared.correlation)) +
    geom_boxplot(aes(fill = Gene)) +
    facet_wrap(~Target_Population, nrow = 2) +
    labs(x = "Comparison population", y = "Spearman correlation squared") +
    theme_bw() + scale_fill_manual(values = c("#CC0000", "#0000FF")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  pdf(paste("GWAS_Plots/An.Overall.Plot.3.Box.plot.GWAS.Correlations.20250509.pdf",sep = ""), width = 12)
    print(p3)
  dev.off()
  

