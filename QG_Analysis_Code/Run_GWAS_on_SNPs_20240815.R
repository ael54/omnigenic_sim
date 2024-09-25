#Run a GS for each subpopulation
list.of.subpopulation.traits <- list(directional.subpopulation.trait.10.pct, directional.subpopulation.trait.20.pct, 
                                     disruptive.subpopulation.trait.10.pct, disruptive.subpopulation.trait.20.pct,
                                     stabilizing.subpopulation.trait.10.pct, stabilizing.subpopulation.trait.20.pct,
                                     directional.subpopulation.trait.10.prev.gen, directional.subpopulation.trait.20.prev.gen,
                                     disruptive.subpopulation.trait.10.prev.gen,disruptive.subpopulation.trait.20.prev.gen,
                                     stabilizing.subpopulation.trait.10.prev.gen,stabilizing.subpopulation.trait.20.prev.gen)

list.of.subpopulation.SNP <- list(directional.subpopulation.10.pct.SNPs, directional.subpopulation.20.pct.SNPs,
                                  disruptive.subpopulation.10.pct.SNPs, disruptive.subpopulation.20.pct.SNPs,
                                  stabilizing.subpopulation.10.pct.SNPs, stabilizing.subpopulation.20.pct.SNPs,
                                  directional.subpopulation.10.prev.gen.SNPs, directional.subpopulation.20.prev.gen.SNPs,
                                  disruptive.subpopulation.10.prev.gen.SNPs, disruptive.subpopulation.20.prev.gen.SNPs,
                                  stabilizing.subpopulation.10.prev.gen.SNPs, stabilizing.subpopulation.20.prev.gen.SNPs)

names.of.subpopulations <- c("Directional.selection.10.pct","Directional.selection.20.pct",
                             "Disruptive.selection.10.pct", "Disruptive.selection.20.pct",
                             "Stabilizing.selection.10.pct","Stabilizing.selection.20.pct",
                             "Directional.selection.10.pct.prev.gen","Directional.selection.20.pct.prev.gen",
                             "Disruptive.selection.10.pct.prev.gen", "Disruptive.selection.20.pct.prev.gen",
                             "Stabilizing.selection.10.pct.prev.gen","Stabilizing.selection.20.pct.prev.gen")


add.effect.estimates.core.SNPs <- list(Directional.selection.10.pct = NA, Directional.selection.20.pct = NA,
                                       Disruptive.selection.10.pct = NA, Disruptive.selection.20.pct = NA,
                                       Stabilizing.selection.10.pct = NA, Stabilizing.selection.20.pct = NA,
                                       Directional.selection.10.pct.prev.gen = NA,Directional.selection.20.pct.prev.gen = NA,
                                       Disruptive.selection.10.pct.prev.gen = NA, Disruptive.selection.20.pct.prev.gen = NA,
                                       Stabilizing.selection.10.pct.prev.gen = NA, Stabilizing.selection.20.pct.prev.gen = NA)

add.effect.estimates.peripheral.SNPs <- list(Directional.selection.10.pct = NA, Directional.selection.20.pct = NA,
                                             Disruptive.selection.10.pct = NA, Disruptive.selection.20.pct = NA,
                                             Stabilizing.selection.10.pct = NA, Stabilizing.selection.20.pct = NA,
                                             Directional.selection.10.pct.prev.gen = NA,Directional.selection.20.pct.prev.gen = NA,
                                             Disruptive.selection.10.pct.prev.gen = NA, Disruptive.selection.20.pct.prev.gen = NA,
                                             Stabilizing.selection.10.pct.prev.gen = NA, Stabilizing.selection.20.pct.prev.gen = NA)



for(eye in 1:length(names.of.subpopulations)){
  this.myY <- data.frame(row.names(list.of.subpopulation.traits[[eye]]),
                         list.of.subpopulation.traits[[eye]])
  this.myGD <- data.frame(row.names(list.of.subpopulation.SNP[[eye]]),list.of.subpopulation.SNP[[eye]]) #(genotypes)
  this.myGM <- data.frame(paste("X", the.physical.map.of.SNPs[,1],sep = ""),
                          the.physical.map.of.SNPs[,2],
                          the.physical.map.of.SNPs[,4])#(map locations)
  this.name.of.GAPIT.directory <- paste("GWAS.of.", names.of.subpopulations[eye],
                                        ".SNPs",sep = "")
  
  
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
  SNP.and.add.EE <- data.frame(myGAPIT$GWAS$SNP,myGAPIT$GWAS$Chromosome,
                               myGAPIT$GWAS$`Position `,
                               myGAPIT$GWAS$Effect.Est)
  
  
  
  #Obtain the row numbers of SNPs that are within 0.05 cM of each QTN
  list.of.row.numbers <- NULL
  for(jay in 1:nrow(four.genetic.values.omni.core.peri.coreperi$core.genes)){
    # Extract the chromosome and bp position of the start site
    row.number.of.QTN <- which(the.physical.map.of.QTLs$id == four.genetic.values.omni.core.peri.coreperi$core.genes$core.genes[jay])
    this.chr.start <- as.numeric(the.physical.map.of.QTLs$chr[row.number.of.QTN])
    this.bp.start <- the.physical.map.of.QTLs$pos[row.number.of.QTN] - 0.05
    
    # Extract the chromosome and bp position of the stop site
    this.chr.stop <- as.numeric(the.physical.map.of.QTLs$chr[row.number.of.QTN])
    this.bp.stop <- the.physical.map.of.QTLs$pos[row.number.of.QTN] + 0.05    
    
    #Identify the row numbers of the desired SNPs
    
    #In most cases the genomic region of interest will be on the same chromosome
    # Therefore the following code will be used
    if(this.chr.start == this.chr.stop){
      these.row.numbers <- which((as.numeric(SNP.and.add.EE$myGAPIT.GWAS.Chromosome)==this.chr.start) 
                                 & (as.numeric(SNP.and.add.EE$myGAPIT.GWAS..Position..) >= this.bp.start) 
                                 & (as.numeric(as.numeric(SNP.and.add.EE$myGAPIT.GWAS.Chromosome)==this.chr.stop) 
                                 & (as.numeric(SNP.and.add.EE$myGAPIT.GWAS..Position..) <= this.bp.stop)))
    }else{
      these.row.numbers.part.1 <- which((as.numeric(SNP.and.add.EE$myGAPIT.GWAS.Chromosome)==this.chr.start)
                                        & (as.numeric(SNP.and.add.EE$myGAPIT.GWAS..Position..) >= this.bp.start)) 
      these.row.numbers.part.2 <- which((as.numeric(as.numeric(SNP.and.add.EE$myGAPIT.GWAS.Chromosome)==this.chr.stop)) 
                                        & (as.numeric(SNP.and.add.EE$myGAPIT.GWAS..Position..) <= this.bp.stop)) 
      these.row.numbers <- c(these.row.numbers.part.1, these.row.numbers.part.2)
    }#end if(this.chr.start == this.chr.stop)
    
    # Append them to a list of positions (eye.e. row numbers)
    list.of.row.numbers <- c(list.of.row.numbers, these.row.numbers)
  }#End for(jay in 1:nrow(the.physical.map.of.QTLs))
  
  
  
  
  
  #Obtain the SNP and allelic effect estimates of the core QTN
  add.EE.of.core.SNPs <-SNP.and.add.EE[list.of.row.numbers,]
  
  add.EE.of.peripheral.SNPs <-SNP.and.add.EE[-list.of.row.numbers,]
  
  add.effect.estimates.core.SNPs[[eye]] <- add.EE.of.core.SNPs
  
  
  add.effect.estimates.peripheral.SNPs[[eye]] <- add.EE.of.peripheral.SNPs
} # End for(eye in 1:length(names.of.subpopulations))

#Calculate the Spearman Rank Correlations of core SNPs.
the.pop.1 <- NULL
the.pop.2 <- NULL
the.spearman.rank <- NULL
for(eye in 1:(length(names.of.subpopulations)-1)){
  for(jay in (eye+1):length(names.of.subpopulations)){
    subpopulation.a <- add.effect.estimates.core.SNPs[[eye]]
    names.of.subpopulation.a <- names.of.subpopulations[eye]
    
    subpopulation.b <- add.effect.estimates.core.SNPs[[jay]]
    names.of.subpopulation.b <- names.of.subpopulations[jay]    
    
    subpopulations.combined <- merge(subpopulation.a, subpopulation.b, 
                                     by.x = "myGAPIT.GWAS.SNP",
                                     by.y = "myGAPIT.GWAS.SNP")
    
    this.spearman.rank.correlation <- cor(subpopulations.combined$myGAPIT.GWAS.Effect.Est.x,
                                          subpopulations.combined$myGAPIT.GWAS.Effect.Est.y,
                                          method = "spearman", use = "complete.obs")
    the.pop.1 <- c(the.pop.1, names.of.subpopulation.a)
    the.pop.2 <- c(the.pop.2, names.of.subpopulation.b)
    the.spearman.rank <- c(the.spearman.rank,  this.spearman.rank.correlation)
    
  }# for(jay in eye:length(names.of.subpopulations))
}#End for(eye in 1:length(names.of.subpopulations))

spearman.correlations.between.core.SNPs <- data.frame(the.pop.1,
                                                      the.pop.2,
                                                      the.spearman.rank)




#Obtain the median correlation across all pairs
median.spearman.correlations.between.core.SNPs<- median(the.spearman.rank)

#Obtain median correlation within each types of selection

keep.this.row <- which(substr(spearman.correlations.between.core.SNPs$the.pop.1,start = 1, stop = 3) ==
                         substr(spearman.correlations.between.core.SNPs$the.pop.2,start = 1, stop = 3))
do.not.keep.this.row <-  unique(c(grep("prev.gen", spearman.correlations.between.core.SNPs$the.pop.1),                       
                                  grep("prev.gen", spearman.correlations.between.core.SNPs$the.pop.2)))
keep.this.row.updated <- keep.this.row[-which(keep.this.row %in% do.not.keep.this.row)]
#same.types.of.selection <- c(same.types.of.selection, keep.this.row)

median.core.spearman.correlation.within.types.of.selection.SNPs  <- 
  median(spearman.correlations.between.core.SNPs$the.spearman.rank[keep.this.row.updated])

#Plot median correlation with previous generation. I am rewriting some code so you can follow
# what I am doing
keep.this.row <- which(substr(spearman.correlations.between.core.SNPs$the.pop.1,start = 1, stop = 3) ==
                         substr(spearman.correlations.between.core.SNPs$the.pop.2,start = 1, stop = 3))
also.keep.this.row.10.pop.1 <- grep(".10.", spearman.correlations.between.core.SNPs$the.pop.1)
also.keep.this.row.10.pop.2 <- grep(".10.", spearman.correlations.between.core.SNPs$the.pop.2)
also.keep.this.row.10 <- also.keep.this.row.10.pop.1[which(also.keep.this.row.10.pop.1 %in% also.keep.this.row.10.pop.2)] 

also.keep.this.row.20.pop.1 <- grep(".20.", spearman.correlations.between.core.SNPs$the.pop.1)
also.keep.this.row.20.pop.2 <- grep(".20.", spearman.correlations.between.core.SNPs$the.pop.2)
also.keep.this.row.20 <- also.keep.this.row.20.pop.1[which(also.keep.this.row.20.pop.1 %in% also.keep.this.row.20.pop.2)] 

also.keep.this.row.10.or.20 <- c(also.keep.this.row.10, also.keep.this.row.20)
keep.this.row.updated.for.current.and.prev.gen <- keep.this.row[which(keep.this.row %in% also.keep.this.row.10.or.20)]

median.core.spearman.correlation.within.prev.gen.SNPs  <- 
  median(spearman.correlations.between.core.SNPs$the.spearman.rank[keep.this.row.updated.for.current.and.prev.gen])




#Obtain median correlation between all subpopulations that are different types of selection
keep.this.row <- which(substr(spearman.correlations.between.core.SNPs$the.pop.1,start = 1, stop = 3) !=
                         substr(spearman.correlations.between.core.SNPs$the.pop.2,start = 1, stop = 3))
do.not.keep.this.row <-  unique(c(grep("prev.gen", spearman.correlations.between.core.SNPs$the.pop.1),                       
                                  grep("prev.gen", spearman.correlations.between.core.SNPs$the.pop.2)))
keep.this.row.updated <- keep.this.row[-which(keep.this.row %in% do.not.keep.this.row)]

median.core.spearman.correlation.different.types.of.selection.SNPs <- 
  median(spearman.correlations.between.core.SNPs$the.spearman.rank[keep.this.row.updated])






#Calculate the Spearman Rank Correlations of Peripheral SNPs.
the.pop.1 <- NULL
the.pop.2 <- NULL
the.spearman.rank <- NULL
for(eye in 1:(length(names.of.subpopulations)-1)){
  for(jay in (eye+1):length(names.of.subpopulations)){
    subpopulation.a <- add.effect.estimates.peripheral.SNPs[[eye]]
    names.of.subpopulation.a <- names.of.subpopulations[eye]
    
    subpopulation.b <- add.effect.estimates.peripheral.SNPs[[jay]]
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

spearman.correlations.between.peripheral.SNPs <- data.frame(the.pop.1,
                                                            the.pop.2,
                                                            the.spearman.rank)

#Obtain the median correlation across all pairs
median.spearman.correlations.between.peripheral.SNPs <- median(the.spearman.rank)

#Obtain median correlation within each types of selection

keep.this.row <- which(substr(spearman.correlations.between.peripheral.SNPs$the.pop.1,start = 1, stop = 3) ==
                         substr(spearman.correlations.between.peripheral.SNPs$the.pop.2,start = 1, stop = 3))
do.not.keep.this.row <-  unique(c(grep("prev.gen", spearman.correlations.between.peripheral.SNPs$the.pop.1),                       
                                  grep("prev.gen", spearman.correlations.between.peripheral.SNPs$the.pop.2)))
keep.this.row.updated <- keep.this.row[-which(keep.this.row %in% do.not.keep.this.row)]
#same.types.of.selection <- c(same.types.of.selection, keep.this.row)

median.peripheral.spearman.correlation.within.types.of.selection.SNPs  <- 
  median(spearman.correlations.between.peripheral.SNPs$the.spearman.rank[keep.this.row.updated])

#Plot median correlation with previous generation. I am rewriting some code so you can follow
# what I am doing
keep.this.row <- which(substr(spearman.correlations.between.peripheral.SNPs$the.pop.1,start = 1, stop = 3) ==
                         substr(spearman.correlations.between.peripheral.SNPs$the.pop.2,start = 1, stop = 3))
also.keep.this.row.10.pop.1 <- grep(".10.", spearman.correlations.between.peripheral.SNPs$the.pop.1)
also.keep.this.row.10.pop.2 <- grep(".10.", spearman.correlations.between.peripheral.SNPs$the.pop.2)
also.keep.this.row.10 <- also.keep.this.row.10.pop.1[which(also.keep.this.row.10.pop.1 %in% also.keep.this.row.10.pop.2)] 

also.keep.this.row.20.pop.1 <- grep(".20.", spearman.correlations.between.peripheral.SNPs$the.pop.1)
also.keep.this.row.20.pop.2 <- grep(".20.", spearman.correlations.between.peripheral.SNPs$the.pop.2)
also.keep.this.row.20 <- also.keep.this.row.20.pop.1[which(also.keep.this.row.20.pop.1 %in% also.keep.this.row.20.pop.2)] 

also.keep.this.row.10.or.20 <- c(also.keep.this.row.10, also.keep.this.row.20)
keep.this.row.updated.for.current.and.prev.gen <- keep.this.row[which(keep.this.row %in% also.keep.this.row.10.or.20)]

median.peripheral.spearman.correlation.within.prev.gen.SNPs  <- 
  median(spearman.correlations.between.peripheral.SNPs$the.spearman.rank[keep.this.row.updated.for.current.and.prev.gen])




#Obtain median correlation between all subpopulations that are different types of selection
keep.this.row <- which(substr(spearman.correlations.between.peripheral.SNPs$the.pop.1,start = 1, stop = 3) !=
                         substr(spearman.correlations.between.peripheral.SNPs$the.pop.2,start = 1, stop = 3))
do.not.keep.this.row <-  unique(c(grep("prev.gen", spearman.correlations.between.peripheral.SNPs$the.pop.1),                       
                                  grep("prev.gen", spearman.correlations.between.peripheral.SNPs$the.pop.2)))
keep.this.row.updated <- keep.this.row[-which(keep.this.row %in% do.not.keep.this.row)]

median.peripheral.spearman.correlation.different.types.of.selection.SNPs <- 
  median(spearman.correlations.between.peripheral.SNPs$the.spearman.rank[keep.this.row.updated])



