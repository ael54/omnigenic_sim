

#Read in phenotypic data

my.Y.validation <- data.frame(row.names(list.of.subpopulation.traits[[this.validation.set]]),
                   list.of.subpopulation.traits[[this.validation.set]])
colnames(my.Y.validation) <- c("taxa", "Trait")

if(length(this.training.set)==1){
  my.Y.training <- data.frame(row.names(list.of.subpopulation.traits[[this.training.set]]),
                                list.of.subpopulation.traits[[this.training.set]])
  colnames(my.Y.training) <- c("taxa", "Trait")
  }else{
  counter <- 0
  for(j in this.training.set){
    this.my.Y.training <- data.frame(row.names(list.of.subpopulation.traits[[j]]),
                                list.of.subpopulation.traits[[j]])
    colnames(this.my.Y.training) <- c("taxa", "Trait")
    if(counter == 0){
      my.Y.training <- this.my.Y.training
    }else{
      my.Y.training <- rbind(my.Y.training, this.my.Y.training)
    }#end if(j == 1)
    counter <- counter+1
  }#end for(j in 1:length(this.training.set))
}#end else

my.Y.for.pipeline <- rbind(my.Y.validation, my.Y.training)
this.trait.name <- colnames(my.Y.for.pipeline)[2]

#Read in the genotypic data
my.G.validation <- data.frame(row.names(list.of.subpopulation.SNPs[[this.validation.set]]),
                              list.of.subpopulation.SNPs[[this.validation.set]])
colnames(my.G.validation)[1] <- "taxa"

if(length(this.training.set)==1){
  my.G.training <- data.frame(row.names(list.of.subpopulation.SNPs[[this.training.set]]),
                              list.of.subpopulation.SNPs[[this.training.set]])
  colnames(my.G.training)[1] <- "taxa"
}else{
  counter <- 0
  for(j in this.training.set){
    this.my.G.training <- data.frame(row.names(list.of.subpopulation.SNPs[[j]]),
                                list.of.subpopulation.SNPs[[j]])
    colnames(this.my.G.training)[1] <- "taxa"
    if(counter == 0){
      my.G.training <- this.my.G.training
    }else{
      my.G.training <- rbind(my.G.training, this.my.G.training)
    }#end if(j == 1)
    counter <- counter + 1
  }#end for(j in 1:length(this.training.set))
}#end else



my.G.for.pipeline <- rbind(my.G.validation, my.G.training)

#########################
#Read in a list of all lines that are in the validation population
taxa.in.validation <- my.Y.validation[,1]

#Temporary code: save my.Y, my.G.whole.genome, and taxa.in.282 into 
# an R workspace to save some time required to read in "large" files


dir.create("Results.GBLUP.Test.20240327")

#Merge the genotypic data to the phenotypic data
#####################Turn this into a function
Y = my.Y.for.pipeline
Geno = my.G.for.pipeline
traitname = this.trait.name
path.for.results = "Results.GBLUP.Test.20240327/"
seed.number = 999


  
 
  y <- as.matrix(Y[,-1])
  
  
  G <- as.matrix(my.G.for.pipeline[,-1], nrow(y), ncol(my.G.for.pipeline[,-1]))
  
  G <- G - 1
  
  #Set up the object CV, which needs to go along for the ride in rrBLUP
  CV=Y[,1:2]
  CV[,2]=1
  colnames(CV)=c("taxa","overall")
  cv <- (as.matrix(CV[,-1]))
  
  taxa.names <- CV[,1]
  
  #Calculate the kinship matrix in rrBLUP
  A1 <- A.mat(G,shrink=TRUE)
  ##############End the function over here - it should ouput: taxa names, A1, and y



#Partition out the merged data so that the 282 is the training set, and the Ames panel is the validation set
pred <- which(taxa.names %in% taxa.in.validation)
  
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








