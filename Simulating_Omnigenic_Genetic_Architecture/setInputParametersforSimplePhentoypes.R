#Initialize the parameters for the function
 
  this.number.of.trait.reps <- 1
  this.broad.sense.H2 <- 0.9
  this.number.of.core.genes <- 7
  this.sigma.2.N <- 1
  this.sigma.2.C <- factor.A*this.sigma.2.N
  this.sigma.2.CC <- factor.B*this.sigma.2.C
  this.sigma.2.NN <- factor.C*this.sigma.2.N 
  this.sigma.2.CN <- factor.D*this.sigma.2.N
  this.seed.number.core.vs.perhiperal <- 2*this.rep
  this.seed.number.core.add <- sample(-1000000:1000000, 1)
  this.seed.number.peripheral.add <- sample(-1000000:1000000, 1)
  this.seed.number.core.core.epi <- sample(-1000000:1000000, 1)
  this.seed.number.peri.peri.epi <- sample(-1000000:1000000, 1)
  this.seed.number.core.peri.epi <- sample(-1000000:1000000, 1)
  this.snps.are.in.columns <- FALSE
  this.output.dir.name <- paste("Factor.A.", factor.A,
                                "Factor.B.", factor.B,
                                "Factor.C.", factor.C,
                                "Factor.D.", factor.D,
                                "Rep.", rep,
                                sep = "")


