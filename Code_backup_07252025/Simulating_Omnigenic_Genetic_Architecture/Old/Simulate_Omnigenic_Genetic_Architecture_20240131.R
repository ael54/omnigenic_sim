#Please see Basic_Idea/Sketch_of_Simulating_Omngenic_Model_20240130..pptx
# for details on how I want to simulate an omnigenic model

#Clear your workspace
rm(list = ls())

#Set the working directory
setwd("/Users/alipka/Library/CloudStorage/Box-Box/Sabbatical_Roslin_Institute/Simulating_Omnigenic_Genetic_Architecture/")
home.dir <- getwd()
#Read in the simplePHENOTYPES library
library(simplePHENOTYPES)
#This is the data set that is included in simplePHENOTYPES. 
# This is the Goodman-Buckler maize diveristy panel (Flint-Garcia et al. 2005)
# and the markers are the from the 55K Illmunia Chip, first published in 
# Cook et al. (2012)

data("SNP55K_maize282_maf04")

####TEMPORARY
# SNP55K_maize282_maf04 has 10,650 SNPs. I am going to randomly\
# select 500 SNPs so that this code can get written (and them sort
# them in the order that they appear), and typos
# can be addressed quickly efficiently.

set.seed(777)
random.sample.of.500.markers <- sample(1:nrow(SNP55K_maize282_maf04), 500)
random.sample.of.500.markers <- random.sample.of.500.markers[order(random.sample.of.500.markers)]

#Input parameters. Note: I am "flipping" the order of this.sigma.2.N
# and this.sigma.2.C (relative to how we intuitively derived this genetic 
# architecture) so that we can make essentially make all 
# variance components proportional to /sigma^2_N

this.input.SNPs <- SNP55K_maize282_maf04[random.sample.of.500.markers, ]
this.broad.sense.H2 <- 0.9
this.number.of.core.genes <- 7
this.sigma.2.N <- 1
this.sigma.2.C <- 2*this.sigma.2.N
this.sigma.2.NN <- 0.05*this.sigma.2.N 
this.sigma.2.CC <- 0.05*this.sigma.2.C
this.sigma.2.CN <- 0.05*this.sigma.2.N
this.seed.number <- 999
this.snps.are.in.columns <- FALSE

#####Values that will eventually go into a function
input.SNPs <- this.input.SNPs
broad.sense.H2 <- this.broad.sense.H2
number.of.core.genes <- this.number.of.core.genes
sigma.2.C <- this.sigma.2.C
sigma.2.N <- this.sigma.2.N
sigma.2.CC <- this.sigma.2.CC
sigma.2.NN <- this.sigma.2.NN
sigma.2.CN <- this.sigma.2.CN 
seed.number <- this.seed.number
snps.are.in.columns <- this.snps.are.in.columns

#Randomly select which SNPs will be the core genes
if(!snps.are.in.columns){
  total.number.of.SNPs <- nrow(input.SNPs)
}else{
  total.number.of.SNPs <- ncol(input.SNPs)
}

#Randomly select core genes 
set.seed(seed.number)
index.of.core.genes <- sample(1:total.number.of.SNPs, number.of.core.genes)

#Assign additive effect sizes to core genes
additive.effects.of.core.genes <- rnorm(n = number.of.core.genes, mean = 0, 
                                        sd = sqrt(sigma.2.C))

#Assign additive effect sizes to peripheral genes
number.of.peripheral.genes <- total.number.of.SNPs - 
  number.of.core.genes
additive.effects.of.peripheral.genes <- rnorm(n = number.of.peripheral.genes, 
                                              mean = 0, 
                                              sd = sqrt(sigma.2.N)) 

#Assign pairwise epistatic effects to all pairs of core genes
number.of.epistatic.effects.within.core <- choose(number.of.core.genes, 2)
epistatic.effects.of.core.genes <- rnorm(n = number.of.epistatic.effects.within.core, 
                                              mean = 0, 
                                              sd = sqrt(sigma.2.CC))

#Assign pairwise epistatic effects to all pairs of peripheral genes
number.of.epistatic.effects.within.peripheral <- choose(number.of.peripheral.genes, 2)
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
# and core-peripheral interacions

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


omnigenic.architecture <-  create_phenotypes(
  geno_obj = input.SNPs,
  add_QTN_num = this.add.QTN.num,
  epi_QTN_num = this.epi.QTN.num,
  h2 = broad.sense.H2,
  add_effect = custom.additive.for.simplePHENOTYPES,
  epi_effect = custom.epistatic.for.simplePHENOTYPES,
  ntraits = 1,
  QTN_list = QTN_list,
  rep = 10,
  vary_QTN = FALSE,
  output_format = "wide",
  architecture = "pleiotropic",
  home_dir = home.dir,
  output_dir = "Initial_Attempt_20230201",
  to_r = T,
  sim_method = "custom",
  seed = 217,
  model = "AE"
)
