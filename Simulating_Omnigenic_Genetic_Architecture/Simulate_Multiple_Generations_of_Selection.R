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
                                          SP.within.function = NULL){
  
    for(i in 1:nGenerations){  
      print(paste("---------------Initiating generation ", i, " of ", 
                  type.of.selection, " selection------------", sep = ""))
      
      #Breeding decision for ith breeding program...select individuals
      if(type.of.selection == "Direct"){
        selected.inds.current.generation = selectInd(current.generation, 
                                        this.nSelect, simParam=SP.within.function)
      }
      if(type.of.selection == "Disruptive"){
        #This code is from https://cran.r-project.org/web/packages/AlphaSimR/AlphaSimR.pdf
        squaredDeviation = function(x, optima = 0) (x - optima)^2
        selected.inds.current.generation = selectInd(current.generation, 
                  this.nSelect, simParam=SP.within.function, 
                  trait = squaredDeviation, selectTop = TRUE)
      }
      if(type.of.selection == "Stabilizing"){
        #This code is from https://cran.r-project.org/web/packages/AlphaSimR/AlphaSimR.pdf
        squaredDeviation = function(x, optima = 0) (x - optima)^2
        selected.inds.current.generation = selectInd(current.generation, 
                           this.nSelect, simParam=SP.within.function, 
                           trait = squaredDeviation, selectTop = FALSE) 
      }
      #Make crosses for the next generation
      the.next.gen = randCross(selected.inds.current.generation, 
                                             this.nCross, simParam=SP.within.function)
      
   
      

      #Get the QTLs
      the.next.genQTL.Genotypes <- pullQtlGeno(pop=the.next.gen, trait = 1, asRaw = FALSE, 
                                                                      simParam = SP.within.function)
      
      
      #Format the QTLs so that they can be read into simplePHENOTYPES
      this.input.QTLs <- get.me.my.SNPs.in.hapmap.format(these.SNPs = the.next.genQTL.Genotypes,
                                                                     this.physical.map = this.QTL.Map)
      
      #Simulate traits with the same QTNs and effect sizes for each generation
      #Get your SNPs
      the.next.genSNPs <- pullSnpGeno(the.next.gen, simParam = SP.within.function)
      
     
      
      #The remaining parameters have been updated already when the founder population was started
      
      this.simulated.trait.within.function <- simulate.omnigenic.architecture(input.SNPs = this.input.QTLs,
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
                                                                     snps.are.in.columns = this.snps.are.in.columns,
                                                                     seed.within.simplePHENOTYPES = this.rep,
                                                                     output.directory.name = this.output.directory.name)
      
      #Put the simulated phentoypes from simplePHNEOTYPES back into AlphaSimR
      the.next.gen@pheno <- as.matrix(this.simulated.trait.within.function$this.simulated.trait[,2])
      
      #Call the next generation the current generation to initiate the next loop
      current.generation <- the.next.gen
    
    } #End for(i in 1:nGenerations)
  
  return(current.generation)

} #end cross.stuff.for.a.whole.bunch.of.generations
