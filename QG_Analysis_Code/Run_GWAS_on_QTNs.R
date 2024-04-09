#Run a GWAS using a GLM for each subpopulation
list.of.subpopulation.traits <- list(directional.subpopulation.trait,
                                     disruptive.subpopulation.trait,
                                     stabilizing.subpopulation.trait)
list.of.subpopulation.QTN <- list(directional.subpopulation.QTNs,
                                  disruptive.subpopulation.QTNs,
                                  stabilizing.subpopulation.QTNs)
#list.of.subpopulation.SNPs <- list(directional.subpopulation.SNPs,
#                                   disruptive.subpopulation.SNPs,
#                                   stabilizing.subpopulation.SNPs)
names.of.subpopulations <- c("Directional.selection", "Disruptive.selection",
                             "Stabilizing.selection")


add.effect.estimates.core.QTN <- list(Directional.selection = NA, 
                                      Disruptive.selection =NA, 
                                      Stabilizing.selection = NA)
add.effect.estimates.peripheral.QTN <- list(Directional.selection = NA, 
                                            Disruptive.selection =NA, 
                                            Stabilizing.selection = NA)
for(eye in 1:length(names.of.subpopulations)){
  this.myY <- data.frame(row.names(list.of.subpopulation.traits[[eye]]),
                         list.of.subpopulation.traits[[eye]])
  this.myGD <- data.frame(row.names(list.of.subpopulation.QTN[[eye]]),list.of.subpopulation.QTN[[eye]]) #(genotypes)
  this.myGM <- data.frame(paste("X", the.physical.map.of.QTLs[,1],sep = ""),
                          the.physical.map.of.QTLs[,2],
                          the.physical.map.of.QTLs[,4])#(map locations)
  this.name.of.GAPIT.directory <- paste("GWAS.of.", names.of.subpopulations[eye],
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
  
  add.effect.estimates.core.QTN[[eye]] <- add.EE.of.core.QTLs
  
  
  add.effect.estimates.peripheral.QTN[[eye]] <- add.EE.of.peripheral.QTLs
} # End for(eye in 1:length(names.of.subpopulations))

#Calculate the Spearman Rank Correlations of core genes.
the.pop.1 <- NULL
the.pop.2 <- NULL
the.spearman.rank <- NULL
for(eye in 1:(length(names.of.subpopulations)-1)){
  for(jay in (eye+1):length(names.of.subpopulations)){
    subpopulation.a <- add.effect.estimates.core.QTN[[eye]]
    names.of.subpopulation.a <- names.of.subpopulations[eye]
    
    subpopulation.b <- add.effect.estimates.core.QTN[[jay]]
    names.of.subpopulation.jay <- names.of.subpopulations[jay]    
    
    subpopulations.combined <- merge(subpopulation.a, subpopulation.b, 
                                     by.x = "myGAPIT.GWAS.SNP",
                                     by.y = "myGAPIT.GWAS.SNP")
    
    this.spearman.rank.correlation <- cor(subpopulations.combined$myGAPIT.GWAS.Effect.Est.x,
                                          subpopulations.combined$myGAPIT.GWAS.Effect.Est.y,
                                          method = "spearman", use = "complete.obs")
    the.pop.1 <- c(the.pop.1, eye)
    the.pop.2 <- c(the.pop.2, jay)
    the.spearman.rank <- c(the.spearman.rank,  this.spearman.rank.correlation)
    
  }# for(jay in eye:length(names.of.subpopulations))
}#End for(eye in 1:length(names.of.subpopulations))

spearman.correlations.between.core.QTNs <- data.frame(the.pop.1,
                                                      the.pop.2,
                                                      the.spearman.rank)



#Calculate the Spearman Rank Correlations of perhiperhal genes.
#Calculate the Spearman Rank Correlations of core genes.
the.pop.1 <- NULL
the.pop.2 <- NULL
the.spearman.rank <- NULL
for(eye in 1:(length(names.of.subpopulations)-1)){
  for(jay in (eye+1):length(names.of.subpopulations)){
    subpopulation.a <- add.effect.estimates.peripheral.QTN[[eye]]
    names.of.subpopulation.a <- names.of.subpopulations[eye]
    
    subpopulation.b <- add.effect.estimates.peripheral.QTN[[jay]]
    names.of.subpopulation.jay <- names.of.subpopulations[jay]    
    
    subpopulations.combined <- merge(subpopulation.a, subpopulation.b, 
                                     by.x = "myGAPIT.GWAS.SNP",
                                     by.y = "myGAPIT.GWAS.SNP")
    
    this.spearman.rank.correlation <- cor(subpopulations.combined$myGAPIT.GWAS.Effect.Est.x,
                                          subpopulations.combined$myGAPIT.GWAS.Effect.Est.y,
                                          method = "spearman",  use = "complete.obs")
    the.pop.1 <- c(the.pop.1, eye)
    the.pop.2 <- c(the.pop.2, jay)
    the.spearman.rank <- c(the.spearman.rank,  this.spearman.rank.correlation )
    
  }# for(jay in eye:length(names.of.subpopulations))
}#End for(eye in 1:length(names.of.subpopulations))

spearman.correlations.between.peripheral.QTNs <- data.frame(the.pop.1,
                                                      the.pop.2,
                                                      the.spearman.rank)


