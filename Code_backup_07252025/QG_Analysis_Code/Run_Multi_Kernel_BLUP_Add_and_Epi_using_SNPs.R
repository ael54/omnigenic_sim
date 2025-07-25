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
  
  #Obtain the row numbers of SNPs that are within 0.05 cM of each QTN
  list.of.col.numbers <- NULL
  for(j in 1:nrow(four.traits.omni.core.peri.coreperi$core.genes)){
    # Extract the chromosome and bp position of the start site
    row.number.of.QTN <- which(the.physical.map.of.QTLs$id == four.traits.omni.core.peri.coreperi$core.genes$core.genes[j])
    this.chr.start <- as.numeric(the.physical.map.of.QTLs$chr[row.number.of.QTN])
    this.bp.start <- the.physical.map.of.QTLs$pos[row.number.of.QTN] - 0.05
    
    # Extract the chromosome and bp position of the stop site
    this.chr.stop <- as.numeric(the.physical.map.of.QTLs$chr[row.number.of.QTN])
    this.bp.stop <- the.physical.map.of.QTLs$pos[row.number.of.QTN] + 0.05    
    
    #Identify the row numbers of the desired SNPs
    
    #In most cases the genomic region of interest will be on the same chromosome
    # Therefore the following code will be used
    if(this.chr.start == this.chr.stop){
      these.col.numbers <- which((as.numeric(the.physical.map.of.SNPs$chr)==this.chr.start) 
                                 & (as.numeric(the.physical.map.of.SNPs$pos) >= this.bp.start) 
                                 & (as.numeric(the.physical.map.of.SNPs$chr)==this.chr.stop) 
                                 & (as.numeric(the.physical.map.of.SNPs$pos) <= this.bp.stop))
    }else{
      these.col.numbers.part.1 <- which((as.numeric(the.physical.map.of.SNPs$chr)==this.chr.start)
                                        & (as.numeric(the.physical.map.of.SNPs$pos) >= this.bp.start)) 
      these.col.numbers.part.2 <- which((as.numeric(as.numeric(the.physical.map.of.SNPs$chr)==this.chr.stop)) 
                                        & (as.numeric(the.physical.map.of.SNPs$pos) <= this.bp.stop)) 
      these.col.numbers <- c(these.row.numbers.part.1, these.row.numbers.part.2)
    }#end if(this.chr.start == this.chr.stop)
    
    #Add 1 to all column numbers because the first column of myG is the taxa namess
    these.col.numbers <- these.col.numbers + 1
    
    # Append them to a list of positions (i.e. row numbers)
    list.of.col.numbers <- c(list.of.col.numbers, these.col.numbers)
  }#End for(j in 1:nrow(the.physical.map.of.QTLs))
  
  
  
  
  my.G.for.pipeline.core <- my.G.for.pipeline[,list.of.col.numbers]
  
  my.G.for.pipeline.peripheral <- my.G.for.pipeline[,-c(1,list.of.col.numbers)]
  
  
  #########################
  #Read in a list of all lines that are in the validation population
  taxa.in.validation <- my.Y.validation[,1]
  
  #Temporary code: save my.Y, my.G.whole.genome, and taxa.in.282 into 
  # an R workspace to save some time required to read in "large" files
  
  
  dir.create("Results.multi.kernel.add.Test.2024032")
  
  #Merge the genotypic data to the phenotypic data
  #####################Turn this into a function
  Y = my.Y.for.pipeline
  traitname = this.trait.name
  path.for.results = "Results.multi.kernel.add.Test.20240327/"
  seed.number = 999
  
  
  
  
  y <- as.matrix(Y[,-1])
  
  
  G.core <- as.matrix(my.G.for.pipeline.core, nrow(y), 
                      ncol(my.G.for.pipeline.core))
  G.core <- G.core - 1
  
  G.peripheral <- as.matrix(my.G.for.pipeline.peripheral, nrow(y), 
                      ncol(my.G.for.pipeline.peripheral))
  G.peripheral <- G.peripheral - 1
  
  
  #Set up the object CV, which needs to go along for the ride in rrBLUP
  CV=Y[,1:2]
  CV[,2]=1
  colnames(CV)=c("taxa","overall")
  cv <- (as.matrix(CV[,-1]))
  
  taxa.names <- CV[,1]

 #####The multi-kernel code differes starting here 
  #Partition G into two groups, depending on whether or not they are in the pathway of interest

  
  #Calculate the partitioned kinship matrices (i.e., kernels)
  A1.core <- A.mat(G.core)
  colnames(A1.core) <- rownames(A1.core) <- c(1:nrow(y))
  
  A1.peripheral <- A.mat(G.peripheral)
  colnames(A1.peripheral) <- rownames(A1.peripheral) <- c(1:nrow(y))  
  
  A1.core.core <- A1.core*A1.core
  
  A1.peripheral.peripheral <- A1.peripheral*A1.peripheral
  
  A1.core.peripheral <- A1.core*A1.peripheral
  
  ############################################################################################
 

  
  sample.size <- nrow(y)
  
 
  #have a "for" loop, start it at 0, and end it at 9
  #I am setting up "k" to denote the nubmer of folds - 1. This is done
  # so that the for loop will work correctly.

  

  pred <- which(taxa.names %in% taxa.in.validation)
  
  yNA <- y
  yNA[pred] <- NA
  
  
  #####data2 is for the double-kernel model. The additional columne, gid2 is needed
  #########so that we can distinguish between the random effects corresponding to 
  ########### each fo the two kernels
  data2 <- data.frame(y=yNA,gid=1:length(y),gid2=1:length(y),gid3=1:length(y),
                      gid4 = 1:length(y), gid5 = 1:length(y), cv = cv)
  the.cv.names <- NULL
  for(j in 1:ncol(cv)) the.cv.names <- c(the.cv.names, paste("CV_",j,sep = ""))
  
  colnames(data2) <- c("y","gid", "gid2","gid3","gid4", "gid5", the.cv.names)
  data2$gid <- as.character(data2$gid)
  data2$gid2 <- as.character(data2$gid2)
  data2$gid3 <- as.character(data2$gid3)
  data2$gid4 <- as.character(data2$gid4)
  data2$gid5 <- as.character(data2$gid5)
  
  
  rownames(A1.core) <- 1:nrow(A1.core)
  rownames(A1.peripheral) <- 1:nrow(A1.peripheral)
  rownames(A1.core.core) <- 1:nrow(A1.core.core)
  rownames(A1.peripheral.peripheral) <- 1:nrow(A1.peripheral.peripheral)
  rownames(A1.core.peripheral) <- 1:nrow(A1.core.peripheral)
  
  #Here is where the magic happens - this is where the two-kernel GS model is fitted
  ans.multiple.kernel <- mmer(y~1, random = ~vsr(gid, Gu = A1.core)+vsr(gid2, Gu = A1.peripheral)
                              +vsr(gid3, Gu = A1.core.core)+vsr(gid4, Gu = A1.peripheral.peripheral)
                              +vsr(gid5, Gu = A1.core.peripheral),
                              data = data2, verbose = FALSE)
  GEBVs <- as.data.frame(ans.multiple.kernel$U$`u:gid`$y+
                           ans.multiple.kernel$U$`u:gid2`$y+
                           ans.multiple.kernel$U$`u:gid3`$y+
                           ans.multiple.kernel$U$`u:gid4`$y+
                           ans.multiple.kernel$U$`u:gid5`$y)
  GEBVs  <- data.frame(as.numeric(rownames(GEBVs)), GEBVs )
  GEBVs <- GEBVs[order(GEBVs[,1]),]  
  
 
  #r.gy <- c(r.gy, cor(ans$g.pred,y[pred]))
  r.gy.add.epi.mult.kern <- cor(GEBVs[pred,2], y[pred]) 
  
  
 #End rrblup.tenfoldCV