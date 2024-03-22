

#Read in phenotypic data
#Work in progress
my.Y.validation <- data.frame(row.names(list.of.subpopulation.traits[[this.validation.set]]),
                   list.of.subpopulation.traits[[this.validation.set]])

my.Y.training <- data.frame(row.names(list.of.subpopulation.traits[[this.training.set]]),
                              list.of.subpopulation.traits[[this.training.set]])



####Old code that need to get updated
#Read in genotypic data
my.G.whole.genome <- read.delim("ames_geno.hmp.txt", head = FALSE) 

#Read in a list of all lines in the Ames panel that are also part of the 282
taxa.in.282 <- read.csv("Accessions_in_282_Only.csv", head = TRUE)

#Temporary code: save my.Y, my.G.whole.genome, and taxa.in.282 into 
# an R workspace to save some time required to read in "large" files

save.image("R.workspace.whole.genome.train.282.predict.Ames.Rdata")

my.Y.for.pipeline <- my.Y[, c(1,2)]
this.trait.name <- colnames(my.Y.for.pipeline)[2]
dir.create("Results.202301224")

#Merge the genotypic data to the phenotypic data
#####################Turn this into a function
Y = my.Y.for.pipeline
Geno = my.G.whole.genome
traitname = this.trait.name
path.for.results = "Results.202301224/"
seed.number = 999





  CV=Y[,1:2]
  CV[,2]=1
  colnames(CV)=c("taxa","overall")
  
  hm=GAPIT.HapMap(G = Geno,SNP.effect="Add",SNP.impute="Major")
  
  
  #####################################
  #Obtain the mafs of all SNPs
  
  #Total number of lines
  ns <- nrow(hm$GD)
  
  #Sum of the allele scores for each SNP
  ss <- apply(hm$GD, 2, sum)
  
  #Combine two situations: one where the allele coded as "2" is major; one where "0" is coded as major.
  maf.matrix <- rbind((.5*ss/ns), (1-(0.5*ss/ns)))
  
  #Copy the minor allele frequencies for all SNPs
  maf <- apply(maf.matrix, 2, min)
  
  #Find out which SNPs have MAF < 0.05
  snps.below.0.05.maf <- which(maf < 0.05)
  
  # Remove these SNPs from hm$GD
  
  hm.GD.without.snps.below.0.05.maf <- hm$GD[,-snps.below.0.05.maf]
  
  ###############################
  
  GK <- cbind(hm$GT, hm.GD.without.snps.below.0.05.maf)
  
  qc=GAPIT.QC(Y = Y, GT = hm$GT, CV = CV, GK = GK)
  
  y <- as.matrix(qc$Y[-1])
  
  G <- as.numeric(qc$GK[,-1])
  
  G <- matrix(G, nrow(y), ncol(qc$GK[,-1]))
  
  G <- G - 1
  
  cv <- (as.matrix(qc$CV[,-1]))
  
  taxa.names <- qc$CV[,1]
  
  #Calculate the kinship matrix in rrBLUP
  A1 <- A.mat(G,shrink=TRUE)
  ##############End the function over here - it should ouput: taxa names, A1, and y



#Partition out the merged data so that the 282 is the training set, and the Ames panel is the validation set
pred <- which(!taxa.names %in% taxa.in.282[,1])
  
#Let GAPIT do this - it can read in the data, match all genotype and phenotype data,
#   and calculate an additive genetic relatedness matrix
yNA <- y
yNA[pred] <- NA

#Train the GS model in the 282
data1 <- data.frame(y=yNA,gid=1:length(y), cv = cv)
the.cv.names <- NULL
for(j in 1:ncol(cv)) the.cv.names <- c(the.cv.names, paste("CV_",j,sep = ""))

colnames(data1) <- c("y","gid", the.cv.names)

rownames(A1) <- 1:nrow(A1)
ans1 <- kin.blup(data1,K=A1,geno="gid",pheno="y", covariate = the.cv.names)


#Evaluate prediction accuracy in the Ames panel
r.gy <-  cor(ans1$g[pred], y[pred])








