#Run a GWAS using a GLM for each subpopulation
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



add.effect.estimates.core.QTN <- list(Directional.selection.10.pct = NA, Directional.selection.20.pct = NA,
                                      Disruptive.selection.10.pct = NA, Disruptive.selection.20.pct = NA,
                                      Stabilizing.selection.10.pct = NA, Stabilizing.selection.20.pct = NA,
                                      Directional.selection.10.pct.prev.gen = NA,Directional.selection.20.pct.prev.gen = NA,
                                      Disruptive.selection.10.pct.prev.gen = NA, Disruptive.selection.20.pct.prev.gen = NA,
                                      Stabilizing.selection.10.pct.prev.gen = NA, Stabilizing.selection.20.pct.prev.gen = NA)

add.effect.estimates.peripheral.QTN <- list(Directional.selection.10.pct = NA, Directional.selection.20.pct = NA,
                                            Disruptive.selection.10.pct = NA, Disruptive.selection.20.pct = NA,
                                            Stabilizing.selection.10.pct = NA, Stabilizing.selection.20.pct = NA,
                                            Directional.selection.10.pct.prev.gen = NA,Directional.selection.20.pct.prev.gen = NA,
                                            Disruptive.selection.10.pct.prev.gen = NA, Disruptive.selection.20.pct.prev.gen = NA,
                                            Stabilizing.selection.10.pct.prev.gen = NA, Stabilizing.selection.20.pct.prev.gen = NA)

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
      paste("X", four.genetic.values.omni.core.peri.coreperi$core.genes$core.genes, sep = "")),]
  
  add.EE.of.peripheral.QTLs <-QTL.and.add.EE[-which(
    QTL.and.add.EE$myGAPIT.GWAS.SNP %in% 
      paste("X", four.genetic.values.omni.core.peri.coreperi$core.genes$core.genes, sep = "")),]
  
  add.effect.estimates.core.QTN[[eye]] <- add.EE.of.core.QTLs
  
  
  add.effect.estimates.peripheral.QTN[[eye]] <- add.EE.of.peripheral.QTLs
} # End for(eye in 1:length(names.of.subpopulations))

#Calculate the Spearman Rank Correlations of core genes.
the.pop.1 <- NULL
the.pop.2 <- NULL
the.spearman.rank <- NULL
for(eye in 1:(length(names.of.subpopulations)-1)){
  for(jay in (eye+1):(length(names.of.subpopulations))){
    subpopulation.a <- add.effect.estimates.core.QTN[[eye]]
    names.of.subpopulation.a <- names.of.subpopulations[eye]
    
    subpopulation.b <- add.effect.estimates.core.QTN[[jay]]
    names.of.subpopulation.b <- names.of.subpopulations[jay]    
    
    subpopulations.combined <- merge(subpopulation.a, subpopulation.b, 
                                     by.x = "myGAPIT.GWAS.SNP",
                                     by.y = "myGAPIT.GWAS.SNP")
    
    this.spearman.rank.correlation <- cor(subpopulations.combined$myGAPIT.GWAS.Effect.Est.x,
                                          subpopulations.combined$myGAPIT.GWAS.Effect.Est.y,
                                          method = "spearman", use = "complete.obs")
    the.pop.1 <- c(the.pop.1, names.of.subpopulation.a)
    the.pop.2 <- c(the.pop.2, names.of.subpopulation.b )
    the.spearman.rank <- c(the.spearman.rank,  this.spearman.rank.correlation)
    
  }# for(jay in eye:length(names.of.subpopulations))
}#End for(eye in 1:length(names.of.subpopulations))

spearman.correlations.between.core.QTNs <- data.frame(the.pop.1,
                                                      the.pop.2,
                                                      the.spearman.rank)


#Obtain the median correlation across all pairs
median.spearman.correlations.between.core.QTN <- median(the.spearman.rank)

#Obtain median correlation within each types of selection

keep.this.row <- which(substr(spearman.correlations.between.core.QTNs$the.pop.1,start = 1, stop = 3) ==
                       substr(spearman.correlations.between.core.QTNs$the.pop.2,start = 1, stop = 3))
do.not.keep.this.row <-  unique(c(grep("prev.gen", spearman.correlations.between.core.QTNs$the.pop.1),                       
                         grep("prev.gen", spearman.correlations.between.core.QTNs$the.pop.2)))
keep.this.row.updated <- keep.this.row[-which(keep.this.row %in% do.not.keep.this.row)]
#same.types.of.selection <- c(same.types.of.selection, keep.this.row)

median.core.spearman.correlation.within.types.of.selection.QTN  <- 
median(spearman.correlations.between.core.QTNs$the.spearman.rank[keep.this.row.updated ])

#Plot median correlation with previous generation. I am rewriting some code so you can follow
# what I am doing
keep.this.row <- which(substr(spearman.correlations.between.core.QTNs$the.pop.1,start = 1, stop = 3) ==
                         substr(spearman.correlations.between.core.QTNs$the.pop.2,start = 1, stop = 3))
also.keep.this.row.10.pop.1 <- grep(".10.", spearman.correlations.between.core.QTNs$the.pop.1)
also.keep.this.row.10.pop.2 <- grep(".10.", spearman.correlations.between.core.QTNs$the.pop.2)
also.keep.this.row.10 <- also.keep.this.row.10.pop.1[which(also.keep.this.row.10.pop.1 %in% also.keep.this.row.10.pop.2)] 

also.keep.this.row.20.pop.1 <- grep(".20.", spearman.correlations.between.core.QTNs$the.pop.1)
also.keep.this.row.20.pop.2 <- grep(".20.", spearman.correlations.between.core.QTNs$the.pop.2)
also.keep.this.row.20 <- also.keep.this.row.20.pop.1[which(also.keep.this.row.20.pop.1 %in% also.keep.this.row.20.pop.2)] 

also.keep.this.row.10.or.20 <- c(also.keep.this.row.10, also.keep.this.row.20)
keep.this.row.updated.for.current.and.prev.gen <- keep.this.row[which(keep.this.row %in% also.keep.this.row.10.or.20)]

median.core.spearman.correlation.within.prev.gen.QTN  <- 
  median(spearman.correlations.between.core.QTNs$the.spearman.rank[keep.this.row.updated.for.current.and.prev.gen])




#Obtain median correlation between all subpopulations that are different types of selection
keep.this.row <- which(substr(spearman.correlations.between.core.QTNs$the.pop.1,start = 1, stop = 3) !=
                         substr(spearman.correlations.between.core.QTNs$the.pop.2,start = 1, stop = 3))
do.not.keep.this.row <-  unique(c(grep("prev.gen", spearman.correlations.between.core.QTNs$the.pop.1),                       
                                  grep("prev.gen", spearman.correlations.between.core.QTNs$the.pop.2)))
keep.this.row.updated <- keep.this.row[-which(keep.this.row %in% do.not.keep.this.row)]

median.core.spearman.correlation.different.types.of.selection.QTN <- 
  median(spearman.correlations.between.core.QTNs$the.spearman.rank[keep.this.row.updated])





#Calculate the Spearman Rank Correlations of perhiperhal genes.
the.pop.1 <- NULL
the.pop.2 <- NULL
the.spearman.rank <- NULL
for(eye in 1:(length(names.of.subpopulations)-1)){
  for(jay in (eye+1):length(names.of.subpopulations)){
    subpopulation.a <- add.effect.estimates.peripheral.QTN[[eye]]
    names.of.subpopulation.a <- names.of.subpopulations[eye]
    
    subpopulation.b <- add.effect.estimates.peripheral.QTN[[jay]]
    names.of.subpopulation.b <- names.of.subpopulations[jay]    
    
    subpopulations.combined <- merge(subpopulation.a, subpopulation.b, 
                                     by.x = "myGAPIT.GWAS.SNP",
                                     by.y = "myGAPIT.GWAS.SNP")
    
    this.spearman.rank.correlation <- cor(subpopulations.combined$myGAPIT.GWAS.Effect.Est.x,
                                          subpopulations.combined$myGAPIT.GWAS.Effect.Est.y,
                                          method = "spearman",  use = "complete.obs")
    the.pop.1 <- c(the.pop.1, names.of.subpopulation.a)
    the.pop.2 <- c(the.pop.2, names.of.subpopulation.b)
    the.spearman.rank <- c(the.spearman.rank,  this.spearman.rank.correlation )
    
  }# for(jay in eye:length(names.of.subpopulations))
}#End for(eye in 1:length(names.of.subpopulations))

spearman.correlations.between.peripheral.QTNs <- data.frame(the.pop.1,
                                                      the.pop.2,
                                                      the.spearman.rank)


#Obtain the median correlation across all pairs
median.spearman.correlations.between.peripheral.QTN <- median(the.spearman.rank)

#Obtain median correlation within each types of selection

keep.this.row <- which(substr(spearman.correlations.between.peripheral.QTNs$the.pop.1,start = 1, stop = 3) ==
                         substr(spearman.correlations.between.peripheral.QTNs$the.pop.2,start = 1, stop = 3))
do.not.keep.this.row <-  unique(c(grep("prev.gen", spearman.correlations.between.peripheral.QTNs$the.pop.1),                       
                                  grep("prev.gen", spearman.correlations.between.peripheral.QTNs$the.pop.2)))
keep.this.row.updated <- keep.this.row[-which(keep.this.row %in% do.not.keep.this.row)]
#same.types.of.selection <- c(same.types.of.selection, keep.this.row)

median.peripheral.spearman.correlation.within.types.of.selection.QTN  <- 
  median(spearman.correlations.between.peripheral.QTNs$the.spearman.rank[keep.this.row.updated])

#Plot median correlation with previous generation. I am rewriting some code so you can follow
# what I am doing
keep.this.row <- which(substr(spearman.correlations.between.peripheral.QTNs$the.pop.1,start = 1, stop = 3) ==
                         substr(spearman.correlations.between.peripheral.QTNs$the.pop.2,start = 1, stop = 3))
also.keep.this.row.10.pop.1 <- grep(".10.", spearman.correlations.between.peripheral.QTNs$the.pop.1)
also.keep.this.row.10.pop.2 <- grep(".10.", spearman.correlations.between.peripheral.QTNs$the.pop.2)
also.keep.this.row.10 <- also.keep.this.row.10.pop.1[which(also.keep.this.row.10.pop.1 %in% also.keep.this.row.10.pop.2)] 

also.keep.this.row.20.pop.1 <- grep(".20.", spearman.correlations.between.peripheral.QTNs$the.pop.1)
also.keep.this.row.20.pop.2 <- grep(".20.", spearman.correlations.between.peripheral.QTNs$the.pop.2)
also.keep.this.row.20 <- also.keep.this.row.20.pop.1[which(also.keep.this.row.20.pop.1 %in% also.keep.this.row.20.pop.2)] 

also.keep.this.row.10.or.20 <- c(also.keep.this.row.10, also.keep.this.row.20)
keep.this.row.updated.for.current.and.prev.gen <- keep.this.row[which(keep.this.row %in% also.keep.this.row.10.or.20)]

median.peripheral.spearman.correlation.within.prev.gen.QTN  <- 
  median(spearman.correlations.between.peripheral.QTNs$the.spearman.rank[keep.this.row.updated.for.current.and.prev.gen])




#Obtain median correlation between all subpopulations that are different types of selection
keep.this.row <- which(substr(spearman.correlations.between.peripheral.QTNs$the.pop.1,start = 1, stop = 3) !=
                         substr(spearman.correlations.between.peripheral.QTNs$the.pop.2,start = 1, stop = 3))
do.not.keep.this.row <-  unique(c(grep("prev.gen", spearman.correlations.between.peripheral.QTNs$the.pop.1),                       
                                  grep("prev.gen", spearman.correlations.between.peripheral.QTNs$the.pop.2)))
keep.this.row.updated <- keep.this.row[-which(keep.this.row %in% do.not.keep.this.row)]

median.peripheral.spearman.correlation.different.types.of.selection.QTN <- 
  median(spearman.correlations.between.peripheral.QTNs$the.spearman.rank[keep.this.row.updated])


