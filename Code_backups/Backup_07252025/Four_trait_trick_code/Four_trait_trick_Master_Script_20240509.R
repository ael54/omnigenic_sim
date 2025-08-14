#Set the working directory

#The first time you run this, please be sure to install
# AlphaSimR and simplePHENOTYPES


#################################################
#Partition out the which loci and effects for each of the four traits
#Trait 1: the whole thing
#Trait 2: only core effects
#Trait 3: only peripheral effects
#Trait 4: only core-by-peripheral epistatic effects



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
                   current.generation = TrainPop,
                   traits = four.traits.omni.core.peri.coreperi$this.simulated.trait,
                   core.genes = four.traits.omni.core.peri.coreperi$core.genes,
                   peripheral.genes = four.traits.omni.core.peri.coreperi$peripheral.genes,
                   core.core.epistasis = four.traits.omni.core.peri.coreperi$core.core.epistasis,
                   peri.peri.epistasis = four.traits.omni.core.peri.coreperi$peri.peri.epistasis,
                   core.peri.epistasis = four.traits.omni.core.peri.coreperi$core.peri.epistasis,
                   calculate.epi.dev = FALSE)


var.covar.of.breeding.values <- popVar(the.breeding.values)
