#Run a GWAS using a GLM for each subpopulation
list.of.subpopulation.traits <- list(directional.subpopulation.trait,
                                     disruptive.subpopulation.trait,
                                     stabilizing.subpopulation.trait)
#list.of.subpopulation.QTN <- list(directional.subpopulation.QTNs,
#                                  disruptive.subpopulation.QTNs,
#                                 stabilizing.subpopulation.QTNs)
list.of.subpopulation.SNPs <- list(directional.subpopulation.SNPs,
                                   disruptive.subpopulation.SNPs,
                                   stabilizing.subpopulation.SNPs)
names.of.subpopulations <- c("Directional.selection", "Disruptive.selection",
                             "Stabilizing.selection")


add.effect.estimates.core.SNPs <- list(Directional.selection = NA, 
                                      Disruptive.selection =NA, 
                                      Stabilizing.selection = NA)
add.effect.estimates.peripheral.SNPs <- list(Directional.selection = NA, 
                                            Disruptive.selection =NA, 
                                            Stabilizing.selection = NA)
for(i in 1:length(names.of.subpopulations)){
  this.myY <- data.frame(row.names(list.of.subpopulation.traits[[i]]),
                         list.of.subpopulation.traits[[i]])
  this.myGD <- data.frame(row.names(list.of.subpopulation.SNPs[[i]]),list.of.subpopulation.SNPs[[i]]) #(genotypes)
  this.myGM <- data.frame(paste("X", the.physical.map.of.SNPs[,1],sep = ""),
                          the.physical.map.of.SNPs[,2],
                          the.physical.map.of.SNPs[,4])#(map locations)
  this.name.of.GAPIT.directory <- paste("GWAS.of.", names.of.subpopulations[i],
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
  for(j in 1:nrow(this.simulated.trait$core.genes)){
    # Extract the chromosome and bp position of the start site
    row.number.of.QTN <- which(the.physical.map.of.QTLs$id == this.simulated.trait$core.genes$core.genes[j])
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
    
    # Append them to a list of positions (i.e. row numbers)
    list.of.row.numbers <- c(list.of.row.numbers, these.row.numbers)
  }#End for(j in 1:nrow(the.physical.map.of.QTLs))
  
  
  
  
  
  #Obtain the SNP and allelic effect estimates of the core QTN
  add.EE.of.core.SNPs <-SNP.and.add.EE[list.of.row.numbers,]
  
  add.EE.of.peripheral.SNPs <-SNP.and.add.EE[-list.of.row.numbers,]
  
  add.effect.estimates.core.SNPs[[i]] <- add.EE.of.core.SNPs
  
  
  add.effect.estimates.peripheral.SNPs[[i]] <- add.EE.of.peripheral.SNPs
} # End for(i in 1:length(names.of.subpopulations))

#Calculate the Spearman Rank Correlations of core SNPs.
the.pop.1 <- NULL
the.pop.2 <- NULL
the.spearman.rank <- NULL
for(i in 1:(length(names.of.subpopulations)-1)){
  for(j in (i+1):length(names.of.subpopulations)){
    subpopulation.a <- add.effect.estimates.core.SNPs[[i]]
    names.of.subpopulation.a <- names.of.subpopulations[i]
    
    subpopulation.b <- add.effect.estimates.core.SNPs[[j]]
    names.of.subpopulation.j <- names.of.subpopulations[j]    
    
    subpopulations.combined <- merge(subpopulation.a, subpopulation.b, 
                                     by.x = "myGAPIT.GWAS.SNP",
                                     by.y = "myGAPIT.GWAS.SNP")
    
    this.spearman.rank.correlation <- cor(subpopulations.combined$myGAPIT.GWAS.Effect.Est.x,
                                          subpopulations.combined$myGAPIT.GWAS.Effect.Est.y,
                                          method = "spearman", use = "complete.obs")
    the.pop.1 <- c(the.pop.1, i)
    the.pop.2 <- c(the.pop.2, j)
    the.spearman.rank <- c(the.spearman.rank,  this.spearman.rank.correlation)
    
  }# for(j in i:length(names.of.subpopulations))
}#End for(i in 1:length(names.of.subpopulations))

spearman.correlations.between.core.SNPs <- data.frame(the.pop.1,
                                                      the.pop.2,
                                                      the.spearman.rank)



#Calculate the Spearman Rank Correlations of Peripheral SNPs.
the.pop.1 <- NULL
the.pop.2 <- NULL
the.spearman.rank <- NULL
for(i in 1:(length(names.of.subpopulations)-1)){
  for(j in (i+1):length(names.of.subpopulations)){
    subpopulation.a <- add.effect.estimates.peripheral.SNPs[[i]]
    names.of.subpopulation.a <- names.of.subpopulations[i]
    
    subpopulation.b <- add.effect.estimates.peripheral.SNPs[[j]]
    names.of.subpopulation.j <- names.of.subpopulations[j]    
    
    subpopulations.combined <- merge(subpopulation.a, subpopulation.b, 
                                     by.x = "myGAPIT.GWAS.SNP",
                                     by.y = "myGAPIT.GWAS.SNP")
    
    this.spearman.rank.correlation <- cor(subpopulations.combined$myGAPIT.GWAS.Effect.Est.x,
                                          subpopulations.combined$myGAPIT.GWAS.Effect.Est.y,
                                          method = "spearman",  use = "complete.obs")
    the.pop.1 <- c(the.pop.1, i)
    the.pop.2 <- c(the.pop.2, j)
    the.spearman.rank <- c(the.spearman.rank,  this.spearman.rank.correlation )
    
  }# for(j in i:length(names.of.subpopulations))
}#End for(i in 1:length(names.of.subpopulations))

spearman.correlations.between.peripheral.SNPs <- data.frame(the.pop.1,
                                                            the.pop.2,
                                                            the.spearman.rank)


