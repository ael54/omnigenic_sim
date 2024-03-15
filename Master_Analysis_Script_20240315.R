#Please see Slide 8: of /Lab_Presentations/Lab_Meeting_20240208/Lipka_A_Omnigenic_20240208.pptx
# for a visualization of what I want to do with this simualtion

#Clear your workspace
rm(list = ls())

#Set the working directory
setwd("/Users/alipka/Library/CloudStorage/Box-Box/Sabbatical_Roslin_Institute/R_workspace/Sabbatical_Project")
home.dir <- getwd()



###Read in the R workspace of simulated data
load("Test_Run_small_data_set_Quick_Haplo_instead_of_runMacs2_20240314.Rdata")


#####Read in all of the packages that are necessary
#Read in prerequiste libaries for GAPIT
library('MASS')
library(multtest)
library(gplots)

#Read in GAPIT
setwd("Scripts_Necessary_for_GAPIT")
source("GAPIT_EMMA source code.txt")
source("GAPIT_Code_from_Internet_20120411_Allelic_Effect.r")
setwd(home.dir)

################################################################################
################################################################################
###GWAS phase
#For the QTLs

#Run a GWAS using a GLM for each subpopulation
list.of.subpopulation.traits <- list(directional.subpopulation.trait,
                                     disruptive.subpopulation.trait,
                                     stabilizing.subpopulation.trait)
list.of.subpopulation.QTN <- list(directional.subpopulation.QTNs,
                                  disruptive.subpopulation.QTNs,
                                  stabilizing.subpopulation.QTNs)
list.of.subpopulation.SNPs <- list(directional.subpopulation.SNPs,
                                   disruptive.subpopulation.SNPs,
                                   stabilizing.subpopulation.SNPs)
names.of.subpopulations <- c("Directional.selection", "Disruptive.selection",
                             "Stabilizing.selection")


add.effect.estimates.core.QTN <- list(Directional.selection = NA, 
                                      Disruptive.selection =NA, 
                                      Stabilizing.selection = NA)
add.effect.estimates.peripheral.QTN <- list(Directional.selection = NA, 
                                            Disruptive.selection =NA, 
                                            Stabilizing.selection = NA)
for(i in 1:length(names.of.subpopulations)){
    this.myY <- data.frame(row.names(list.of.subpopulation.traits[[i]]),
                           list.of.subpopulation.traits[[i]])
    this.myGD <- data.frame(row.names(list.of.subpopulation.QTN[[i]]),list.of.subpopulation.QTN[[i]]) #(genotypes)
    this.myGM <- data.frame(paste("X", the.physical.map.of.QTLs[,1],sep = ""),
                            the.physical.map.of.QTLs[,2],
                            the.physical.map.of.QTLs[,4])#(map locations)
    this.name.of.GAPIT.directory <- paste("GWAS.of.", names.of.subpopulations[i],
                                          ".QTNs",sep = "")


    dir.create(paste(this.name.of.GAPIT.directory))
    setwd(paste(this.name.of.GAPIT.directory))
    myGAPIT <- GAPIT(
      Y=this.myY,				#This is phenotype data
      GD=this.myGD,				#This is genotype data,set it to NULL with multiple genotype files
      GM=this.myGM,
      PCA.total = 0,
      group.from=1,		#Was 232	#Lower bound for number of group
      group.to=1,			#Upper bound for number of group
      group.by=1,				#range between 1 and number of individuals, smaller the finer
      SNP.impute = "Major",
      SNP.MAF = 0.00,
      cutOff = 0.00,
    )
    setwd(home.dir)
    
    #Obtain the SNP and allelic effect estimates from GAPIT
    QTL.and.add.EE <- data.frame(myGAPIT$GWAS$SNP,
                                     myGAPIT$GWAS$Effect.Est)
    
    #Obtain the SNP and allelic effect estimates of the core QTN
    add.EE.of.core.QTLs <-QTL.and.add.EE[which(
        QTL.and.add.EE$myGAPIT.GWAS.SNP %in% 
         paste("X", this.simulated.trait$core.genes$core.genes, sep = "")),]
    
    add.EE.of.peripheral.QTLs <-QTL.and.add.EE[-which(
      QTL.and.add.EE$myGAPIT.GWAS.SNP %in% 
        paste("X", this.simulated.trait$core.genes$core.genes, sep = "")),]
    
    add.effect.estimates.core.QTN[[i]] <- add.EE.of.core.QTLs

    
    add.effect.estimates.peripheral.QTN[[i]] <- add.EE.of.peripheral.QTLs
} # End for(i in 1:length(names.of.subpopulations))

#Calculate the Spearman Rank Correlations of core genes.
the.pop.1 <- NULL
the.pop.2 <- NULL
the.spearman.rank <- NULL
for(i in 1:(length(names.of.subpopulations)-1)){
  for(j in (i+1):length(names.of.subpopulations)){
    subpopulation.a <- add.effect.estimates.core.QTN[[i]]
    names.of.subpopulation.a <- names.of.subpopulations[i]

    subpopulation.b <- add.effect.estimates.core.QTN[[j]]
    names.of.subpopulation.j <- names.of.subpopulations[j]    
    
    subpopulations.combined <- merge(subpopulation.a, subpopulation.b, 
                                    by.x = "myGAPIT.GWAS.SNP",
                                    by.y = "myGAPIT.GWAS.SNP")
    
    this.spearman.rank.correlation <- cor(subpopulations.combined$myGAPIT.GWAS.Effect.Est.x,
                                          subpopulations.combined$myGAPIT.GWAS.Effect.Est.y,
                                          method = "spearman")
    the.pop.1 <- c(the.pop.1, i)
    the.pop.2 <- c(the.pop.2, j)
    the.spearman.rank <- c(the.spearman.rank,  this.spearman.rank.correlation )
    
  }# for(j in i:length(names.of.subpopulations))
}#End for(i in 1:length(names.of.subpopulations))

spearman.correlations.between.core.QTNs <- data.frame(the.pop.1,
                                                        the.pop.2,
                                                        the.spearman.rank)




    #Obtain the SNP and allelic effect estimates of the peripheral QTN
#
  #Quantify the Spearman rank correlation between allelic effects
  #  of core QTLs between subpopulations


  #Quantify the Spearman rank correlation between allelic effects
  #  of peripheral QTLs between subpopulations


#For the SNPs
  
  #Run a GWAS using a GLM for each subpopulation
  #  
  
  #Quantify the Spearman rank correlation between allelic effects
  #  of SNPs in the vicinity of core QTLs between subpopulations
  
  
  #Quantify the Spearman rank correlation between allelic effects
  #  of SNPs in the vicinity of peripheral QTLs between subpopulations


################################################################################
################################################################################
###Estimating variance components phase

#For the QTLs




#For the SNPs





################################################################################
################################################################################
###Genomic prediction across subpopulations phase

#For the QTLs




#For the SNPs




