library(AlphaSimR)

# Create haplotypes for an initial population
founderPop = runMacs(nInd=1000, nChr=1, segSites=100)

# Add a trait with additive and additive-by-additive effects
SP = SimParam$
  new(founderPop)$
  addTraitAE(100, relAA=0.3, var=100, mean=10)

# Create the initial population
pop = newPop(founderPop)


### Calculate contribution of additive effects to variance components
# This is really simple because additive effects only contribution to breeding values
# Also, the average effect of an allele substitution equals the additive effect

# Extract QTL genotypes
M = pullQtlGeno(pop)

# Center genotypes for calculation of breeding value, this sucks off the mean
M_centered = scale(M, scale=FALSE)

# Get additive (average) effects
a = SP$traits[[1]]@addEff

# Calculate breeding values
BV = M_centered %*% a 

# Calculating the contribution to the trait mean from the additive effects and the intercept
mu = mean((M-1)%*%a) + SP$traits[[1]]@intercept

### Calculate contribution of additive-by-additive epistatic effects to variance components
# They will contribute to both the breeding values and the additive-by-additive epistatic deviations

# Get additive-by-additive effects
aa = SP$traits[[1]]@epiEff

# Calculate frequencies for all genotype classes
# This is needed for calculating the average effect of an allele substitution
freq = matrix(0, 
              nrow=3, # Three genotype classes for diploids (0,1,2)
              ncol=ncol(M))

for(i in seq_len(3)){
  freq[i,] = colSums(M == (i-1)) / nrow(M)
}

# Create a variable to hold the additive-by-additive epistatic deviations
# Making it a matrix to match BV above
AA = matrix(0, nrow=nrow(BV))

# Loop over all pairs of interacting loci
#
# I'll calculate the average effect for each locus to find breeding values
# I'll then calculate the additive-by-additive epistatic deviations using the 
# lack-of-fit after accounting for the breeding values 
# This approach works because there is no other form of epistasis in the simulation.
for(i in seq_len(nrow(aa))){
  
  # Extract positions for loci and the epistatic effect
  loc1 = aa[i,1] # Locus 1
  loc2 = aa[i,2] # Locus 2
  eff = aa[i,3] # Effect
  
  # Construct centered genotypes for both loci
  x = c(-1,0,1)
  x1 = x - sum(freq[,loc1]*x)
  x2 = x - sum(freq[,loc2]*x)
  
  # Calculate marginal means for locus 1 and center
  # Values correspond to genotypes 0, 1, and 2
  y1 = c((freq[1,loc2]-freq[3,loc2])*eff,
         0,
         (freq[3,loc2]-freq[1,loc2])*eff)
  y1 = y1 - sum(y1*freq[,loc1])
  
  # Calculate marginal means for locus 2 and center
  # Values correspond to genotypes 0, 1, and 2
  y2 = c((freq[1,loc1]-freq[3,loc1])*eff,
         0,
         (freq[3,loc1]-freq[1,loc1])*eff)
  y2 = y2 - sum(y2*freq[,loc2])
  
  # Calculate average effects for both loci
  # Average effects are just regression coefficients
  alpha = c(
    sum(freq[,loc1]*y1*x1) / sum(freq[,loc1]*x1^2),
    sum(freq[,loc2]*y2*x2) / sum(freq[,loc2]*x2^2)
  )
  
  # Account for a potential divide by zero case
  if(any(is.na(alpha))) alpha[is.na(alpha)] = 0
  
  # Calculate breeding values
  tmpBv = M_centered[, c(loc1, loc2)] %*% alpha
  
  # Add breeding values to the global breeding values
  # At this point the values could be kept separate to partition variance due 
  # to different types of QTL
  BV = BV + tmpBv
  
  # Calculate epistatic deviations using lack-of-fit
  # This calculation could be skipped if you don't care about epistatic variance
  # It takes considerably more time than the previous calculations
  
  # Construct a 3x3 table for the epistatic deviations
  # The rows will represent the 0, 1, and 2 genotypes for locus 1
  # The columns will represent the 0, 1, and 2 genotypes for locus 2
  # This will be done in three parts
  # 1. Construct initial matrix using genetic values
  # 2. Subtract average genetic value using observed genotype frequencies (includes LD)
  # 3. Subtract breeding values
  devTab = tcrossprod(x) * eff
  freqTab = matrix(0, nrow=3, ncol=3)
  for(j in seq_len(nrow(M))){
    freqTab[M[j,loc1]+1, M[j,loc2]+1] = freqTab[M[j,loc1]+1, M[j,loc2]+1] + 1
  }
  freqTab = freqTab / nrow(M)
  tmpMu = sum(freqTab*devTab)
  devTab = devTab - tmpMu
  devTab = sweep(devTab, 1, x1*alpha[1], "-")
  devTab = sweep(devTab, 2, x2*alpha[2], "-")
  
  # Determine which element of the table to access for each individual
  # This is being done to avoid using a double loop
  # To make sense of this trick, note that R saves data in a column major order
  # Thus, a value of 3 would be column 1, row 3 and a value of 4 would be column 2, row 1
  ind = M[,loc1] + 3*M[,loc2] + 1
  
  # Get each individual's epistatic deviation
  tmpAa = devTab[ind]
  
  # Add deviations to global value
  AA = AA + tmpAa
  
  # Add contribution to trait mean
  mu = mu + tmpMu
}

# Calculate using AlphaSimR's built-in functionality for comparison
asr = genParam(pop) 

# Examine mean breeding value (zero by definition), note some numeric instability
mean(BV)
mean(asr$bv)

# Examine mean epistatic deviations (zero by definition), note some numeric instability
mean(AA)
mean(asr$aa)

# Calculate additive genetic variance (variance of breeding values)
sum(BV^2) / length(BV)
asr$varA

# Calculate additive-by-additive genetic variance (variance of epistatic deviations)
sum(AA^2) / length(AA)
asr$varAA





