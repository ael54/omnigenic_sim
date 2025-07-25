#Simulated trait object from the founder population
this.simulated.trait


#Simulated trait object from one generation of directional selection
this.simulated.trait.direct


#Check for equality between core, peripheral, and corresponding epistatic sets
setequal(this.simulated.trait$core.genes, this.simulated.trait.direct$core.genes)
setequal(this.simulated.trait$peripheral.genes, this.simulated.trait.direct$peripheral.genes)
setequal(this.simulated.trait$core.core.epistasis, this.simulated.trait.direct$core.core.epistasis)
setequal(this.simulated.trait$peri.peri.epistasis, this.simulated.trait.direct$peri.peri.epistasis)
setequal(this.simulated.trait$core.peri.epistasis, this.simulated.trait.direct$core.peri.epistasis)
