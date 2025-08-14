#Simulate a trait. We are not actually going to use this trait for anything
# once we have a founder population, so I will hard code in the settings
# from Bancic et al. (2023)
SP$addTraitAG(nQtlPerChr = this.nQtl,
              mean       = 1,
              var        = 1,
              varEnv     = 1,
              varGxE     = 4)

# Collect pedigree
SP$setTrackPed(TRUE)


# Create founder parents
Parents = newPop(founderPop)

# Add phenotype reflecting evaluation in EYT
Parents = setPheno(Parents, varE = this.varE.during.burnin, reps = this.repEYT)