##############################################
##############################################
##############################################  
#Created by Alex Lipka on February 21, 2024. This
# Uses functions and objects available in AlphaSimR and simplePHENTOYPES

cross.stuff.for.a.whole.bunch.of.generations <- function(current.generation = NULL,
                                          type.of.selection = NULL,
                                          nSelect = NULL,
                                          nCross = NULL,
                                          nGenerations = NULL,
                                          nQtl = NULL,
                                          SP.within.function = NULL,
                                          trait.var.covar = NULL,
                                          genetic.value.var.covar = NULL,
                                          breeding.value.var.covar = NULL,
                                          trait.values = NULL,
                                          genetic.values = NULL,
                                          breeding.values = NULL,
                                          initial.count.value = 1){
    count.local <- initial.count.value+1
    if(!((length(trait.var.covar) == length(genetic.value.var.covar))&
         (length(trait.var.covar) == length(breeding.value.var.covar)))){
      print("####################The lengths of trait.var.covar, genetic.value.var.covar, and breeding.value.var.covar are not equal. Please fix this.#############################")
  
    }
    for(i in 1:nGenerations){  
      print(paste("---------------Initiating generation ", i, " of ", 
                  type.of.selection, " selection------------", sep = ""))
      
      #Store the current generation as the "previous generation". This
      # will allow us to use info of the previous generation in GWAS and GS
      previous.generation = current.generation
     
       #Breeding decision for ith breeding program...select individuals
      if(type.of.selection == "Direct"){
        selected.inds.current.generation = selectInd(current.generation, 
                                        nSelect, simParam=SP.within.function)
      }
      if(type.of.selection == "Disruptive"){
        #This code is from https://cran.r-project.org/web/packages/AlphaSimR/AlphaSimR.pdf
        squaredDeviation = function(x, optima = 0) (x - optima)^2
        selected.inds.current.generation = selectInd(current.generation, 
                  nSelect, simParam=SP.within.function, 
                  trait = squaredDeviation, selectTop = TRUE)
      }
      if(type.of.selection == "Stabilizing"){
        #This code is from https://cran.r-project.org/web/packages/AlphaSimR/AlphaSimR.pdf
        squaredDeviation = function(x, optima = 0) (x - optima)^2
        selected.inds.current.generation = selectInd(current.generation, 
                           nSelect, simParam=SP.within.function, 
                           trait = squaredDeviation, selectTop = FALSE) 
      }
      #Make crosses for the next generation
      the.next.gen = randCross(selected.inds.current.generation, 
                                             this.nCross, simParam=SP.within.function)
      
   
      

      #Get the QTLs
      the.next.genQTL.Genotypes <- pullQtlGeno(pop=the.next.gen, trait = 1, asRaw = FALSE, 
                                                                      simParam = SP.within.function)
      
      
      #Format the QTLs so that they can be read into simplePHENOTYPES
      this.input.SNPs <- get.me.my.SNPs.in.hapmap.format(these.SNPs = the.next.genQTL.Genotypes,
                                                                     this.physical.map = this.QTL.Map)
      
      #Simulate traits with the same QTNs and effect sizes for each generation
      #Get your SNPs
      the.next.genSNPs <- pullSnpGeno(the.next.gen, simParam = SP.within.function)
      
     
      
      #The remaining parameters have been updated already when the founder population was started
      
      #Simulate the four traits in simplePHENOTYPES
      four.traits.omni.core.peri.coreperi <- simulate.omni.four.trait.trick (input.SNPs = this.input.SNPs,
                                                                             number.of.trait.reps = this.number.of.trait.reps,
                                                                             broad.sense.H2 = this.broad.sense.H2,
                                                                             number.of.core.genes = this.number.of.core.genes,
                                                                             sigma.2.C = this.sigma.2.C,
                                                                             sigma.2.N = this.sigma.2.N,
                                                                             sigma.2.CC = this.sigma.2.CC,
                                                                             sigma.2.NN = this.sigma.2.NN,
                                                                             sigma.2.CN = this.sigma.2.CN, 
                                                                             seed.number.core.vs.perhiperal = this.seed.number.core.vs.perhiperal,
                                                                             seed.number.core.add = this.seed.number.core.add,
                                                                             seed.number.peripheral.add = this.seed.number.peripheral.add,
                                                                             seed.number.core.core.epi = this.seed.number.core.core.epi,
                                                                             seed.number.peri.peri.epi = this.seed.number.peri.peri.epi,
                                                                             seed.number.core.peri.epi = this.seed.number.core.peri.epi,
                                                                             seed.within.simplePHENOTYPES = this.rep,
                                                                             snps.are.in.columns = this.snps.are.in.columns,
                                                                             output.directory.name = this.output.directory.name)
      
      #Simulate the genetic values of the four traits in simplePHENOTYPES
      four.genetic.values.omni.core.peri.coreperi <- simulate.omni.four.trait.trick (input.SNPs = this.input.SNPs,
                                                                                     number.of.trait.reps = this.number.of.trait.reps,
                                                                                     broad.sense.H2 = 1,
                                                                                     number.of.core.genes = this.number.of.core.genes,
                                                                                     sigma.2.C = this.sigma.2.C,
                                                                                     sigma.2.N = this.sigma.2.N,
                                                                                     sigma.2.CC = this.sigma.2.CC,
                                                                                     sigma.2.NN = this.sigma.2.NN,
                                                                                     sigma.2.CN = this.sigma.2.CN, 
                                                                                     seed.number.core.vs.perhiperal = this.seed.number.core.vs.perhiperal,
                                                                                     seed.number.core.add = this.seed.number.core.add,
                                                                                     seed.number.peripheral.add = this.seed.number.peripheral.add,
                                                                                     seed.number.core.core.epi = this.seed.number.core.core.epi,
                                                                                     seed.number.peri.peri.epi = this.seed.number.peri.peri.epi,
                                                                                     seed.number.core.peri.epi = this.seed.number.core.peri.epi,
                                                                                     seed.within.simplePHENOTYPES = this.rep,
                                                                                     snps.are.in.columns = this.snps.are.in.columns,
                                                                                     output.directory.name = this.output.directory.name)
      
      
      #Use popVar to calculate the variance-covariance matrix of the traits
      the.trait.values <- as.matrix(four.traits.omni.core.peri.coreperi$this.simulated.trait[,-c(1,6)])
      var.covar.of.trait.values <- popVar(the.trait.values)
     
      
      #Use popVar to calculate the variance-covariance matrix of the genetic values
      the.genetic.values <- as.matrix(four.genetic.values.omni.core.peri.coreperi$this.simulated.trait[,-c(1,6)])
      var.covar.of.genetic.values <- popVar(the.genetic.values)

      
      #Obtain the breeding values
      the.breeding.values <- get.me.my.breeding.values(SP.within.function = SP,
                                                       current.generation = the.next.gen,
                                                       traits = four.traits.omni.core.peri.coreperi$this.simulated.trait,
                                                       core.genes = four.traits.omni.core.peri.coreperi$core.genes,
                                                       peripheral.genes = four.traits.omni.core.peri.coreperi$peripheral.genes,
                                                       core.core.epistasis = four.traits.omni.core.peri.coreperi$core.core.epistasis,
                                                       peri.peri.epistasis = four.traits.omni.core.peri.coreperi$peri.peri.epistasis,
                                                       core.peri.epistasis = four.traits.omni.core.peri.coreperi$core.peri.epistasis,
                                                       calculate.epi.dev = FALSE)
      
      
      var.covar.of.breeding.values <- popVar(the.breeding.values)
          
      #Put the simulated phenotypes from simplePHNEOTYPES back into AlphaSimR
      the.next.gen@pheno <- as.matrix(four.traits.omni.core.peri.coreperi$this.simulated.trait[,2])
      
      
      #Add the population variance-covariance matrix of trait.values, genetic.values, and breeding values to the list
      trait.var.covar[[count.local]] <- var.covar.of.trait.values
      names(trait.var.covar)[count.local] <- paste(type.of.selection, ".nSel.",nSelect,".Gen.",i,sep = "") 
      
      genetic.value.var.covar[[count.local]] <- var.covar.of.genetic.values
      names(genetic.value.var.covar)[count.local] <-  paste(type.of.selection, ".nSel.",nSelect, ".Gen.",i,sep = "")
      
      breeding.value.var.covar[[count.local]] <- var.covar.of.breeding.values
      names(breeding.value.var.covar)[count.local] <-  paste(type.of.selection, ".nSel.",nSelect, ".Gen.",i,sep = "")
      
      
      #Add the observed trait.values, genetic.values, and breeding values to the list 
      trait.values[[count.local]] <- the.trait.values 
      names(trait.values)[count.local] <- paste(type.of.selection, ".nSel.",nSelect,".Gen.",i,sep = "") 
      
      genetic.values[[count.local]] <- the.genetic.values
      names(genetic.values )[count.local] <- paste(type.of.selection, ".nSel.",nSelect,".Gen.",i,sep = "") 
      
      breeding.values[[count.local]] <- the.breeding.values
      names(breeding.values)[count.local] <- paste(type.of.selection, ".nSel.",nSelect,".Gen.",i,sep = "") 
      
      #Call the next generation the current generation to initiate the next loop
      current.generation = the.next.gen
      
      
      
      count.local <- count.local+1
    } #End for(i in 1:nGenerations)
  
    #NOTE: trait.values, genetic.values, and breeding.values are a list of such values for all ten generations,
    ## while four.trait.values, four.genetic.values, and four.breeding values are just for the last generation.
  return(list(current.generation = current.generation, 
              previous.generation = previous.generation,
              trait.var.covar = trait.var.covar, 
              genetic.value.var.covar = genetic.value.var.covar,
              breeding.value.var.covar = breeding.value.var.covar,
              trait.values = trait.values,
              genetic.values = genetic.values,
              breeding.values = breeding.values,
              four.trait.values = the.trait.values,
              four.genetic.values = the.genetic.values,
              four.breeding.values = the.breeding.values))

} #end cross.stuff.for.a.whole.bunch.of.generations
