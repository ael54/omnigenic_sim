#Please see Basic_Idea/Sketch_of_Simulating_Omngenic_Model_20240130..pptx
# for details on how I want to simulate an omnigenic model
simulate.omnigenic.architecture <- function(input.SNPs = NULL,
                                            number.of.trait.reps = NULL,
                                            broad.sense.H2 = NULL,
                                            number.of.core.genes = NULL,
                                            sigma.2.C = NULL,
                                            sigma.2.N = NULL,
                                            sigma.2.CC = NULL,
                                            sigma.2.NN = NULL,
                                            sigma.2.CN = NULL, 
                                            seed.number.core.vs.perhiperal = NULL,
                                            seed.number.core.add = NULL,
                                            seed.number.peripheral.add = NULL,
                                            seed.number.core.core.epi = NULL,
                                            seed.number.peri.peri.epi = NULL,
                                            seed.number.core.peri.epi = NULL,
                                            seed.within.simplePHENOTYPES = NULL,
                                            snps.are.in.columns = TRUE,
                                            output.directory.name = NULL){
  #Randomly select which SNPs will be the core genes
  if(!snps.are.in.columns){
    total.number.of.SNPs <- nrow(input.SNPs)
  }else{
    total.number.of.SNPs <- ncol(input.SNPs)
  }
  
  #Randomly select core genes 
  set.seed(seed.number.core.vs.perhiperal)
  index.of.core.genes <- sample(1:total.number.of.SNPs, number.of.core.genes)
  
  #Assign additive effect sizes to core genes
  set.seed(seed.number.core.add)
  additive.effects.of.core.genes <- rnorm(n = number.of.core.genes, mean = 0, 
                                          sd = sqrt(sigma.2.C))
  
  #Assign additive effect sizes to peripheral genes
  number.of.peripheral.genes <- total.number.of.SNPs - 
    number.of.core.genes
  set.seed(seed.number.peripheral.add)
  additive.effects.of.peripheral.genes <- rnorm(n = number.of.peripheral.genes, 
                                                mean = 0, 
                                                sd = sqrt(sigma.2.N)) 
  
  #Assign pairwise epistatic effects to all pairs of core genes
  number.of.epistatic.effects.within.core <- choose(number.of.core.genes, 2)
  set.seed(seed.number.core.core.epi)
  epistatic.effects.of.core.genes <- rnorm(n = number.of.epistatic.effects.within.core, 
                                                mean = 0, 
                                                sd = sqrt(sigma.2.CC))
  
  #Assign pairwise epistatic effects to all pairs of peripheral genes
  number.of.epistatic.effects.within.peripheral <- choose(number.of.peripheral.genes, 2)
  set.seed(seed.number.peri.peri.epi)
  epistatic.effects.of.peripheral.genes <- rnorm(n = number.of.epistatic.effects.within.peripheral, 
                                           mean = 0, 
                                           sd = sqrt(sigma.2.NN))
  
  #Assign pairwise epistatic effects to all pairs of core and peripheral genes
  #####The following number is derived from the numerator of the 
  # probability mass function of a hypergeometric random variable
  ### E.g. suppose we have k core genes and (k-p) peripheral genes. How many possible ways
  ######## can we obtain 1 core gene and 1 periperhal gene
  number.of.epistatic.effects.between.core.and.peripheral <- 
    choose(number.of.core.genes,1)*choose(number.of.peripheral.genes, 1)
    
  set.seed(seed.number.core.peri.epi)
  epistatic.effects.between.core.and.perhpheral <- rnorm(n = number.of.epistatic.effects.between.core.and.peripheral, 
                                           mean = 0, 
                                           sd = sqrt(sigma.2.CN))
  
  
  #Now do the simulations on simplePHENOTYPES
  
  #Put all of the additive and epistatic effects into a single vector
  # of additive effects, followed by a single vector of epistastic effects
  custom.additive.for.simplePHENOTYPES <- list(trait_1 = c(additive.effects.of.core.genes, 
                                            additive.effects.of.peripheral.genes))
  custom.epistatic.for.simplePHENOTYPES <- list(trait_1 = c(epistatic.effects.of.core.genes,
                                             epistatic.effects.of.peripheral.genes,
                                             epistatic.effects.between.core.and.perhpheral))
  
  #specify which SNPs are core genes, and which are periperhal genes
  if(!snps.are.in.columns){
    SNP.names <- input.SNPs$snp
  }else{
    #Do something else to get a list of 
    # markers; this will be figured out soon
  }
  
  core.genes <- SNP.names[index.of.core.genes]
  peripheral.genes <- SNP.names[-index.of.core.genes]
  
  #exhaustively list all possible combinations of two core genes
  core.pairwise.epistasis <- NULL
  for(i in 1:(length(core.genes)-1)){
    for(j in (i+1): length(core.genes)){
      core.pairwise.epistasis <- c(core.pairwise.epistasis, core.genes[i],
                                   core.genes[j])
    }#end for(j in 2: number.of.epistatic.effects.within.core)
  }#end for(i in 1:number.of.epistatic.effects.within.core)
  
  #exhaustively list all possible combinations of two peripheral genes
  peripheral.pairwise.epistasis <- NULL
  for(i in 1:(length(peripheral.genes)-1)){
    for(j in (i+1): length(peripheral.genes)){
      peripheral.pairwise.epistasis <- c(peripheral.pairwise.epistasis, 
                                         peripheral.genes[i],
                                         peripheral.genes[j])
    }#end for(j in 2: number.of.epistatic.effects.within.core)
  }#end for(i in 1:number.of.epistatic.effects.within.core)
  
  
  #exhaustively list all possible combinations of one core and one peripheral genes
  core.peripheral.pairwise.epistasis <- NULL
  for(i in 1:(length(core.genes))){
    for(j in 1 : length(peripheral.genes)){
      core.peripheral.pairwise.epistasis <- c(core.peripheral.pairwise.epistasis, 
                                         core.genes[i],
                                         peripheral.genes[j])
    }#end for(j in 2: number.of.epistatic.effects.within.core)
  }#end for(i in 1:number.of.epistatic.effects.within.core)
  
  ####And now, we have the necessary input parameters to specify which
  # SNPs are core genes, which are peripheral genes, and which pairs
  # of SNPs are core-core interaction, peripheral-peripheral interactions,
  # and core-peripheral interactions
  
  QTN_list <- list()
  QTN_list$add[[1]] <- c(core.genes, peripheral.genes)
  QTN_list$epi[[1]] <- c(core.pairwise.epistasis, peripheral.pairwise.epistasis,
                    core.peripheral.pairwise.epistasis)
  
  
  ####Tying in some loose ends:
  this.add.QTN.num <- length(custom.additive.for.simplePHENOTYPES$trait_1)
  this.epi.QTN.num <- length(custom.epistatic.for.simplePHENOTYPES$trait_1)
  #####This code below based on  
  #https://cran.r-project.org/web/packages/simplePHENOTYPES/vignettes/simplePHENOTYPES.html
  # I am going to adjust this to simulated the genetic architecture that 
  # we want.
  
  print("Yee haw!!! We are simulating some traits now!")
  
  omnigenic.architecture <-  create_phenotypes(
    geno_obj = input.SNPs,
    add_QTN_num = this.add.QTN.num,
    epi_QTN_num = this.epi.QTN.num,
    h2 = broad.sense.H2,
    add_effect = custom.additive.for.simplePHENOTYPES,
    epi_effect = custom.epistatic.for.simplePHENOTYPES,
    ntraits = 1,
    QTN_list = QTN_list,
    rep = number.of.trait.reps,
    vary_QTN = FALSE,
    output_format = "wide",
    architecture = "pleiotropic",
    home_dir = home.dir,
    output_dir = output.directory.name,
    to_r = T,
    sim_method = "custom",
    seed = seed.within.simplePHENOTYPES,
    model = "AE"
  )
  
  #Compile information on which QTNs have which additive and
  # epistatic effects
  core.genes <- data.frame(core.genes =  core.genes, 
                           add.effects = additive.effects.of.core.genes)
  
  peripheral.genes <- data.frame(peripheral.genes = peripheral.genes,
                                 add.effects = additive.effects.of.peripheral.genes)
  
  epistatic.core.effects.for.export <- c(sapply(epistatic.effects.of.core.genes, c,
                                         NA))

  
  core.core.epistasis <- data.frame(core.core.epi = core.pairwise.epistasis,
                                    epi.effects =   epistatic.core.effects.for.export)
  
  
  epistatic.peri.effects.for.export <- c(sapply(epistatic.effects.of.peripheral.genes, c,
                                                NA))
  
  
  peri.peri.epistasis <- data.frame(peri.peri.epi = peripheral.pairwise.epistasis,
                                    epi.effects = epistatic.peri.effects.for.export)
  
  
  epistatic.core.peri.effects.for.export <- c(sapply(epistatic.effects.between.core.and.perhpheral, c,
                                                NA))
  
  core.peri.epistasis <- data.frame(core.peri.epi = core.peripheral.pairwise.epistasis,
                                    epi.effects = epistatic.core.peri.effects.for.export)
  
  
  
  return(list(this.simulated.trait = omnigenic.architecture,
              core.genes = core.genes,
              peripheral.genes = peripheral.genes,
              core.core.epistasis = core.core.epistasis,
              peri.peri.epistasis = peri.peri.epistasis,
              core.peri.epistasis = core.peri.epistasis))
}# End simulate.omnigenic.architecture